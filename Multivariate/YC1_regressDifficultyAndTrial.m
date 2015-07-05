function YC1_regressDifficultyAndTrial(subjs,saveDir)
% function YC1_regressDifficultyAndTrial(subjs)
%
% Option inputs: cells array of subject strings, defaults to all subjects
%                save directory
%
% For each encoding trial, at all timepoints and frequencies, perform a
% regression of the following type:
%
%      power = B0 + B1(difficulty) + B2(learningTrialNumber) + B3(trialNum)
%
% where:
%      
%      difficulty is the average performace factor score across all trials
%      in the entire dataset within a 5 VR radius of the object location
%
%      learningTrialNumber indicates either the first or second learning
%      trial
%
%      trialNum is trial number within the session
%
% This function saves the following files to saveDir/subj
%
%       <subj>_regStatistics.mat
%       <subj>_residuals.mat

    

% save directory
if ~exist('saveDir','var') || isempty(saveDir)
    saveDir = '/data10/scratch/jfm2/YC1/multi/power/regress';
end
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% do bipolar
bipol = 1;

if ~exist('params','var')
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

% see if this was submitted with an open pool
poolobj = gcp('nocreate'); 
if ~isempty(poolobj)
    tic
    parfor s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        runRegress_subj(subjs{s},bipol,params,allErrors,allObjectLocs,saveDir);
    end
    toc    
else
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        runRegress_subj(subjs{s},bipol,params,allErrors,allObjectLocs,saveDir);
        
    end
end
end



function runRegress_subj(subj,bipol,params,allErrors,allObjectLocs,saveDir)

% load tal structure
try
    tal = getBipolarSubjElecs(subj,bipol,1);
catch
    fprintf('Could not load electrode locations for %s. Aborting.\n',subj)
    return
end
    
if ~isfield(tal,'locTag') || ~any(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
    fprintf('No MTL electrodes for %s.\n',subj)
    return
end
nElecs = length(tal);

subjDir = fullfile(saveDir,subj);
if ~exist(subjDir,'dir')
    mkdir(subjDir)
end

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

% get beta1 (avg difficulty). This is based on average error of all test
% trials within a certain number of VR units from the object location,
% excluding the current trial. CURRENTLY USING 5 VR-UNIT RADIUS
avgDiff = NaN(sum(eventsToUse),1);
objLocs = vertcat(events(eventsToUse).objLocs);
for trial = 1:length(avgDiff)
    x = objLocs(trial,1);
    y = objLocs(trial,2);
    near = sqrt((allObjectLocs(:,1) - x).^2 + (allObjectLocs(:,2) - y).^2) < 5;    
    current = ismember(allObjectLocs,objLocs(trial,:),'rows');
    avgDiff(trial) = mean(allErrors(near & ~current));
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
x = (x-repmat(mean(x),size(x,1),1)) ./ repmat(std(x),size(x,1),1);

% loop over electrode
for e = 1:size(powerData,2)
    fprintf('%s: Electrode %d of %d.\n',subj,e,size(powerData,2))
    
    % Init beta matrices, times x freqs
    [beta1,beta2,beta3] = deal(NaN(size(powerData,3),size(powerData,4)));
    
    % Init tstat matrices, times x freqs
    [tstat1,tstat2,tstat3] = deal(NaN(size(powerData,3),size(powerData,4)));
    
    % Init pval matrices, times x freqs
    [pval1,pval2,pval3] = deal(NaN(size(powerData,3),size(powerData,4)));
    
    % Init residuals matrix, events times x freqs.
    resid = NaN(size(powerData,1),size(powerData,3),size(powerData,4));
    
    % time
    for t = 1:size(powerData,3)
        
        % frequency
        for f = 1:size(powerData,4)
            
            % observations = power for each event for the electrode at this
            % time and frequency
            y = powerData(:,e,t,f);
            
            % run regression
            s = regstats(y,x,'linear',{'beta','yhat','r','mse','rsquare','tstat'});
            
            % save output beta
            beta1(t,f) = s.beta(2);
            beta2(t,f) = s.beta(3);
            beta3(t,f) = s.beta(4);
            
            % and tstats
            tstat1(t,f) = s.tstat.t(2);
            tstat2(t,f) = s.tstat.t(3);
            tstat3(t,f) = s.tstat.t(4);   
            
            % and pvals
            pval1(t,f) = s.tstat.pval(2);
            pval2(t,f) = s.tstat.pval(3);
            pval3(t,f) = s.tstat.pval(4);  
            
            % and finally the residuals, which will serve as corrected
            % power values for later analyses
            resid(:,t,f) = s.r;
            
        end % frequency
    end % time
    
    % save all the betas,tstats,pvals to one file
    vars = {'beta1','beta2','beta3','tstat1','tstat2','tstat3','pval1','pval2','pval3','tal'};
    save(fullfile(subjDir,sprintf('%s_elec_%s_regStatistics.mat',subj,tal(e).eNames)),vars{:});
    
    % save the residuals to another file
    resid = single(resid); %#ok<NASGU>
    save(fullfile(subjDir,sprintf('%s_elec_%s_residuals.mat',subj,tal(e).eNames)),'resid');
    
end % electrode


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
    subjStd = distOut.stdBasePow;
    subjPow = distOut.pow;
    
    % a couple sessions weren't zscored, so do it here. I should double
    % check that this is right
    if isempty(subjStd)
        fprintf('Power not zscored for %s. Doing it now.\n',subj)
        
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