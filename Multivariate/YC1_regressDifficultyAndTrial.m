function [ output_args ] = YC1_regressDifficultyAndTrial( params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    

% save directory
saveDir = '/data10/scratch/jfm2/YC1/multi/power/regress';
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% do bipolar
bipol = 1;

if ~exist('params','var') || isempty(params)
    params = multiParams();
    params.timeBins = [];
    params.freqBins = [];
end

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% load all errors
[allErrors,allObjectLocs] = YC1_loadAllSubjErrors;

if exist('pool','var')
    matlabpool(pool,length(subjs)+1)
    tic
    parfor s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        runRegress_subj(subjs{s},bipol,params,allErrors,allObjectLocs,saveDir);
    end
    toc
    matlabpool close
    toc
else
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        runRegress_subj(subjs{s},bipol,params,allErrors,allObjectLocs,saveDir);
        
    end
end
end



function runRegress_subj(subj,bipol,params,allErrors,allObjectLocs,saveDir)

% 
% fname = fullfile(saveDir,[subj '_lasso.mat']);
% if exist(fname,'file')
%     return
% end

% load tal structure
tal = getBipolarSubjElecs(subj,bipol,1,1);
if ~isfield(tal,'locTag') || ~any(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
    fprintf('No MTL electrodes for %s.\n',subj)
    return
end
nElecs = length(tal);

% load events so we can filter into our conditions of interest
config = RAM_config('RAM_YC1');

% Setting time bins for convenience:
tEnds = (config.distributedParams.timeWin:...
    config.distributedParams.timeStep:...
    config.postMS+config.priorMS)-config.priorMS;
tStarts = tEnds - config.distributedParams.timeWin + 1;
config.distributedParams.timeBins = [tStarts' tEnds'];

% load events
[events] = RAM_loadEvents(subj,[],'RAM_YC1','events', config);

% add the test error to the learning trials
events = addErrorField(events);


% convert to the session to an double from a string
session = NaN(1,length(events));
for e = 1:length(session)
    session(e) = str2double(events(e).session);
end

% filter to events of interest
eventsToUse = params.eventFilter(events);
if sum(eventsToUse) < 10
    fprintf('Not enough events for %s.\n',subj)
    return
end

% get time and freq settings
freqBins = params.freqBins;
timeBins = params.timeBins;

% load power for all electrodes
powerData = loadAllPower(tal,subj,events,freqBins,timeBins,config,eventsToUse);


% reorder to be events x electrodes x time x freq.
powerData = permute(powerData,[3 4 2 1]);

% for every electrode, will perform a separate regression at every time
% point and frequency. Predictors will be:
%
%     beta1 = average difficulty of the target location
%     beta2 = learning trial 1 or learning trial 2
%     beta3 = overall trial number, 1...n, n = trial number within session

% get beta1 (avg difficulty). This is based on all average error of all
% test trials within a certain number of VR units from the object location.
avgDiff = NaN(sum(eventsToUse),1);
objLocs = vertcat(events(eventsToUse).objLocs);
for trial = 1:length(avgDiff)
    x = objLocs(trial,1);
    y = objLocs(trial,2);
    near = sqrt((allObjectLocs(:,1) - x).^2 + (allObjectLocs(:,2) - y).^2) < 5;
    avgDiff(trial) = mean(allErrors(near));
end
    
% get beta2 (learning trial 1 or 2 for each test trial)
[~,IA,~] = unique(objLocs,'rows','stable');
learningNum = ones(sum(eventsToUse),1);
learningNum(IA+1) = 2;

% get beta3 (overall trial number within session)
trialNumber = [events(eventsToUse).itemno]'+1;

% create data matrix
x = [avgDiff learningNum trialNumber];

% normalize?
x =(x-repmat(mean(x),size(x,1),1)) ./ repmat(std(x),size(x,1),1);


% loop over electrode
for e = 1:size(powerData,2)
   
    % time 
    for t = 1:size(powerData,3)
       
        for f = 1:size(powerData,4)
            
            % observations = power for each event for the electrode at this
            % time and frequency
            y = powerData(:,e,t,f);
            keyboard
            
        end % frequency    
    end % time
end % electrode
keyboard
end

function powerData = loadAllPower(tal,subj,events,freqBins,timeBins,config,eventsToUse)

doFreqs = 0;
if size(freqBins,1) > 0
    nFreqs = size(freqBins,1);
    doFreqs = 1;
else
    nFreqs = length(config.distributedParams.freQ);
end

doTimes = 0;
if size(timeBins,1) > 0
    nTimes = size(timeBins,1);
    doTimes = 1;
else
    nTimes = size(config.distributedParams.timeBins,1);
end

nEvents = sum(eventsToUse);
nElecs = length(tal);
powerData = NaN(nFreqs,nTimes,nEvents,nElecs);

for e = 1:nElecs
    elecNum = tal(e).channel;
    
    % load power for all sessions. Power should aleady have been
    % created or else error
    [distOut] = RAM_dist_func(subj,[],elecNum,'RAM_YC1','events', 0, ...
        @doNothing, config.distributedFunctionLabel, config.distributedParams, 1, 1, events);
    
    sessInds = distOut.sessInds;
    subjMean = distOut.meanBasePow;
    subjStd = distOut.stdBasePow;
    subjPow = distOut.pow;
    
    % a couple sessions weren't zscored, so do it here. I should double
    % check that this is right
    if isempty(subjStd)
        fprintf('power not zscored for %s\n',subj)
        
        sI = unique(sessInds);
        zpow = NaN(size(subjPow));
        for s = 1:length(sI)
            inds = sessInds == sI(s);
            subjMean = nanmean(squeeze(nanmean(subjPow(:,:,inds),2)),2);
            subjMean = repmat(subjMean,[1 size(subjPow,2), size(subjPow,3)]);
            subjStd = nanstd(squeeze(nanmean(subjPow(:,:,inds),2)),[],2);
            subjStd = repmat(subjStd,[1 size(subjPow,2), size(subjPow,3)]);
            zpow(:,:,inds) = (subjPow(:,:,inds) - subjMean).*subjStd;
        end
        subjPow = zpow;                
    end
    subjPow = subjPow(:,:,eventsToUse);
    
    % average frequencies
    if doFreqs
        tmpPower = NaN(nFreqs,size(subjPow,2),size(subjPow,3));
        for f = 1:nFreqs
            fInds = config.distributedParams.freQ >= freqBins(f,1) & config.distributedParams.freQ < freqBins(f,2);
            tmpPower(f,:,:) = nanmean(subjPow(fInds,:,:),1);
        end
        subjPow = tmpPower;
    end        
    
    % average times
    if doTimes
        tmpPower = NaN(nFreqs,nTimes,size(subjPow,3));
        for t = 1:nTimes
            tInds = config.distributedParams.timeBins(:,1) >= timeBins(t,1) & config.distributedParams.timeBins(:,2) < timeBins(t,2);
            tmpPower(:,t,:) = nanmean(subjPow(:,tInds,:),2);
        end  
        powerData(:,:,:,e) = tmpPower;        
    else
        powerData(:,:,:,e) = subjPow;
    end
end
end

function events = addErrorField(events)
% add testError field
% add inner field (1 = inner region, 0 = outer region)

testInd = strcmp({events.type},'NAV_TEST');
recEvents = events(testInd);
[events.testError] = deal(NaN);
[events.recalled] = deal(NaN);
[events.inner] = deal(NaN);
sessVec = [events.session];
trialVec = [events.blocknum];
for rec = 1:length(recEvents);
  session = recEvents(rec).session;  
  trial = recEvents(rec).blocknum;
  err = recEvents(rec).respPerformanceFactor;
  ind = sessVec == session & trialVec == trial;
  [events(ind).testError] = deal(err);
  [events(ind).inner] = deal(abs(recEvents(rec).objLocs(1)) < 568/30 && abs(recEvents(rec).objLocs(2)) < 7);
end
end