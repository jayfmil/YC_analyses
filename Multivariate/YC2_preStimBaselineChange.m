function [res] = YC2_preStimBaselineChange(subj,params,yc1Data,chanceData,saveDir)
% function [] = YC2_applyWeights(subj,params,saveDir)
%
% Inputs:
%
%       subj - subject string
%     params - params structure
%    savedir - string of path to save directory
%
%
%
% Saves results to

res = [];
try
    
    % make sure params.modelEachTime is true or that the number of timebins
    % is only one, otherwise we can't apply to model to just the one
    % post-stim time bin we are intereseted in
    if params.modelEachTime == 0 && size(params.timeBins,2) > 1
        fprintf('Cannot apply YC1 model...too many time bins for %s.\n',subj);
    end
    
    % Use same electrodes as YC1
    yc1Tal = yc1Data.tal;
    
    % also load tal structure for all subject elecs so we can ensure we can
    % look up the anatomical location
    talAll = getBipolarSubjElecs(subj,1);
    
    % load power parameters
    %%%% WE ARE PROBABLY GOING TO NEED SEPERATE YC1 AND YC2 params?
    powParams = load(fullfile(params.powerPath,'params.mat'));
    
    % Setting time bins for convenience
    tEnds     = (powParams.params.pow.timeWin:powParams.params.pow.timeStep:powParams.params.eeg.durationMS)+powParams.params.eeg.offsetMS;
    tStarts   = tEnds - powParams.params.pow.timeWin+1;
    powParams.timeBins = [tStarts' tEnds'];
    
    % load events
    events = get_sub_events('RAM_YC2',subj);
    
    % add the test error to the learning trials
    events  = addErrorField(events);
    
    % this is stupid
    if strcmp(subj,'R1047D')
        [events.stimLeads] = deal('LOTD3-LOTD4');
    end    
    
    % The level of analysis is unique stim locations
    stimLeads   = unique({events.stimLeads});
    res = [];
    for stimLoc = 1:length(stimLeads)        
        
        res(stimLoc).stimTagName = stimLeads{stimLoc};
        res(stimLoc).stimAnat    = talAll(strcmp({talAll.tagName},res(stimLoc).stimTagName)).locTag;
        fprintf('Processing %s: stim region: %s, %s (%d of %d).\n',subj,res(stimLoc).stimAnat,res(stimLoc).stimTagName,stimLoc,length(stimLeads))
        
        % filter to events with the stim locations
        eventsToUse = strcmp({events.stimLeads},stimLeads{stimLoc});
        
        % further filter to just no practice encoding
%         eventsToUse = eventsToUse & params.eventFilter(events);
        eventsToUse = eventsToUse & (strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD'));%params.eventFilter(events);
        
        % load post stimulation and matched non-stim power for these events
        
        powerData = loadAllPower(yc1Tal,subj,events,params.freqBins,[-1000 -100],powParams,eventsToUse,params);
%         powerData = loadAllPower(yc1Tal,subj,events,params.freqBins,[5100 7000],powParams,eventsToUse,params);
        powerData = permute(powerData,[3 1 2 4]);
        
        % indices of stim and non-stim trials
        stimEvents    = [events(eventsToUse).isStim]==1;
        nonStimEvents = [events(eventsToUse).isStim]==0;
        
        [resTmp] = computePostStimStats(powerData,yc1Data,chanceData,stimEvents,nonStimEvents,events,eventsToUse);
        res(stimLoc).deltaRR         = resTmp.deltaRR;
        res(stimLoc).deltaRR_bin     = resTmp.deltaRR_binary;
        res(stimLoc).yc1Score        = resTmp.yc1Score;        
        res(stimLoc).numNZweights    = resTmp.numNZweights;
        res(stimLoc).nonStimBaselineChange = resTmp.nonStimBaselineChange;
        res(stimLoc).stimBaselineChange = resTmp.stimBaselineChange;
                        
        % Do 1000 permutattion of stim non-stim labels
        stimLabels           = [events(eventsToUse).isStim];
        deltaRR_perm         = NaN(1000,3);
        deltaRR_bin_perm     = NaN(1000,3);

        
        for i = 1:1000              
            fprintf('Processing %s: stim region: %s, %s (%d of %d). Permuation %d of 1000.\n',subj,res(stimLoc).stimAnat,res(stimLoc).stimTagName,stimLoc,length(stimLeads),i)
            Ytmp           = reshape(stimLabels,2,[]);
            randOrder      = randperm(length(Ytmp));
            stimLabelsRand = reshape(Ytmp(:,randOrder),[],1);
            stimEvents     = stimLabelsRand == 1;
            nonStimEvents  = stimLabelsRand == 0;
            
            % this is inefficient, don't need to redo glmval each time...
            resTmp = computePostStimStats(powerData,yc1Data,chanceData,stimEvents,nonStimEvents,events,eventsToUse);              
            
            deltaRR_perm(i,:)         = resTmp.deltaRR;
            deltaRR_bin_perm(i,:)     = resTmp.deltaRR_binary;          
        end
       
        res(stimLoc).deltaRR_perm         = deltaRR_perm;
        res(stimLoc).deltaRR_bin_perm     = deltaRR_bin_perm;
       
        
        
    end            

    fname = fullfile(saveDir,[subj '_YC2_preStimChange.mat']);
    save(fname,'res');
catch e
    fname = fullfile(saveDir,[subj '_YC2_preStimChange_error.mat']);
    save(fname,'e')
end

function [res] = computePostStimStats(powerData,yc1Data,chanceData,stimEvents,nonStimEvents,events,eventsToUse)
% following Youssef's terminology, compute recall relative to
% non-stim baseline. He uses percent recall, here I use performance
% factor change
res = [];

testErrors  = 1-[events(eventsToUse).testError];

% also make binary
recalled = testErrors>median(testErrors);


% size of feature matrix
nElecs = size(powerData,4);
nFreqs = size(powerData,2);

% Choose most significnat timebin from YC1. Instead of looping
% over all times? If we have a tie, use most sig with highest AUC
p = mean(repmat(yc1Data.AUC,[size(chanceData.auc_all,1), 1]) > chanceData.auc_all);
maxP = max(p);
AUC = yc1Data.AUC;
AUC(p~=maxP) = NaN;
[~,timeToUse] = max(AUC);
res.yc1Score = maxP;

% ALWAYS USE ENC PERIOD. MAKE THIS AN OPTION
timeToUse    = strcmpi(yc1Data.params.timeBinLabels,'Enc');
res.yc1Score = p(timeToUse);

% average weights across all across all YC1 training folds
A                = mean(horzcat(yc1Data.res(timeToUse).A{:}),2);
intercept        = mean([yc1Data.res(timeToUse).intercept{:}]);
res.numNZweights = nnz(A);

% predict YC1 time bin weights applied to YC2 time bin
B1 = [intercept;A];

% predict with YC1 time bin weights applied to YC2 post period
X = reshape(squeeze(powerData(:,:,1,:)),size(powerData,1),nFreqs*nElecs);
preProb = glmval(B1,X,'logit');

[preProbSort,inds] = sort(preProb);
stimEventsSort = stimEvents(inds);
testErrorsSort = testErrors(inds);
recSort = recalled(inds);

baselinePerf = mean(recalled(nonStimEvents));

% now bin the sorted recall vector
start = 1:(length(preProbSort)/3):length(preProbSort);
stop = [start(2:end)-1 length(preProbSort)];
deltaRR = NaN(1,length(stop));
deltaRR_binary = NaN(1,length(stop));

nonStimBaselineChange = NaN(1,length(stop));
stimBaselineChange = NaN(1,length(stop));

for r = 1:length(stop)
    
    probBin = preProbSort(start(r):stop(r));
    stimBin = stimEventsSort(start(r):stop(r));
    testErrorsBin = testErrorsSort(start(r):stop(r));
    recBin = recSort(start(r):stop(r));
    
    preProbStimBin = probBin(stimBin);
    preProbNonStimBin = probBin(~stimBin);
    
    testErrorsStimBin = testErrorsBin(stimBin);
    testErrorsNonStimBin = testErrorsBin(~stimBin);
    
    recStimBin = recBin(stimBin);
    recNonStimBin = recBin(~stimBin);    
    
    deltaRR(r) = 100*((nanmean(testErrorsStimBin) - nanmean(testErrorsNonStimBin))/nanmean(testErrorsNonStimBin));
    deltaRR_binary(r) = 100*((nanmean(recStimBin) - nanmean(recNonStimBin))/nanmean(recNonStimBin));   
        
    nonStimBaselineChange(r) = 100*((nanmean(recNonStimBin)-baselinePerf)/baselinePerf);   
    stimBaselineChange(r) = 100*((nanmean(recStimBin)-baselinePerf)/baselinePerf);
        
end
res.deltaRR = deltaRR;
res.deltaRR_binary = deltaRR_binary;
res.nonStimBaselineChange =nonStimBaselineChange;
res.stimBaselineChange = stimBaselineChange;



function powerData = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params)

nFreqs = size(freqBins,1);
nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
nElecs = length(tal);
powerData = NaN(nFreqs,nTimes,nEvents,nElecs);

% when loading power, use either original power or power with effect of
% trial number removed.
powField = 'pow';
if params.useCorrectedPower
    powField = 'powCorr';
end

for e = 1:nElecs
    elecNum = tal(e).channel;
    
    basePath  = params.powerPath;%'/data10/scratch/jfm2/RAM/biomarker/power/';
    subjPath  = fullfile(basePath,subj);
    sessions = unique([events.session]);
    subjPow  = [];
    for s = 1:length(sessions)
        fname = fullfile(subjPath,'RAM_YC2_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'_pre.mat']);
        sessPow = load(fname);
        subjPow = cat(3,subjPow,sessPow.sessOutput.(powField));
    end
    
    if length(eventsToUse) ~= size(subjPow,3)
        fprintf('Number of events does not match size of power matrix for %s!.\n',subj)
        return
    end
    subjPow = subjPow(:,:,eventsToUse);
    
    % average frequencies
    if nFreqs ~= length(powParams.params.pow.freqs)
      tmpPower = NaN(nFreqs,size(subjPow,2),size(subjPow,3));
      for f = 1:nFreqs
        fInds = powParams.params.pow.freqs >= freqBins(f,1) & powParams.params.pow.freqs < freqBins(f,2);
        tmpPower(f,:,:) = nanmean(subjPow(fInds,:,:),1);
      end
      subjPow = tmpPower;
    end
    
    % average times
%     tmpPower = NaN(nFreqs,nTimes,size(subjPow,3));
%     for t = 1:nTimes
%         tInds = powParams.timeBins(:,1) >= timeBins(t,1) & powParams.timeBins(:,2) < timeBins(t,2);
%         tmpPower(:,t,:) = nanmean(subjPow(:,tInds,:),2);
%     end
    powerData(:,:,:,e) = subjPow;
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




