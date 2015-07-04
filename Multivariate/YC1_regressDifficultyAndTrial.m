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

if exist('pool','var')
    matlabpool(pool,length(subjs)+1)
    tic
    parfor s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        runRegress_subj(subjs{s},bipol,params,saveDir);
    end
    toc
    matlabpool close
    toc
else
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        runRegress_subj(subjs{s},bipol,params,saveDir);
        
    end
end
end



function runRegress_subj(subj,bipol,params,saveDir)

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
powerData = permute(powerData,[3 1 2 4]);

end

function powerData = loadAllPower(tal,subj,events,freqBins,timeBins,config,eventsToUse)

if size(freqBins,1) > 0
    nFreqs = size(freqBins,1);
else
    nFreqs = length(config.distributedParams.freQ);
end
if size(timeBins,1) > 0
    nTimes = size(timeBins,1);
else
    nTimes = 4
end
keyboard
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
    if nFreqs > 0
        tmpPower = NaN(nFreqs,size(subjPow,2),size(subjPow,3));
        for f = 1:nFreqs
            fInds = config.distributedParams.freQ >= freqBins(f,1) & config.distributedParams.freQ < freqBins(f,2);
            tmpPower(f,:,:) = nanmean(subjPow(fInds,:,:),1);
        end
        subjPow = tmpPower;
    end        
    
    % average times
    if nTimes > 0
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