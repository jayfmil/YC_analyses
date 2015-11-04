function [res] = YC2_postStimBaselineChange(subj,params,yc1Data,chanceData,stimToUse,timeToUse,saveDir)
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
    powParams = load(fullfile(params.powerPath,'params.mat'));
    
    % Setting time bins for convenience
    tEnds     = (powParams.params.pow.timeWin:powParams.params.pow.timeStep:powParams.params.eeg.durationMS)+powParams.params.eeg.offsetMS;
    tStarts   = tEnds - powParams.params.pow.timeWin+1;
    powParams.timeBins = [tStarts' tEnds'];
    
    % load yc1 events
    eventsYC1 = get_sub_events('RAM_YC1',subj);
    yc1thresh = median([eventsYC1(strcmp({eventsYC1.type},'NAV_TEST')).respPerformanceFactor]);    
    
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
        session = [events.session];
        eventsToUse = eventsToUse & (strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD'));%params.eventFilter(events);
        
        % and finally filter to stim of interest (first, second, or both)
        if any(strcmp({'first','second'},stimToUse))
            [~,firstIdx,~] = unique([session' [events.blocknum]'],'rows','first');
            if strcmp(stimToUse,'first')
                eventsToUse(firstIdx+1) = false;
            elseif strcmp(stimToUse,'second')
                eventsToUse(firstIdx) = false;
            end
        end
                       
        % load post stimulation and matched non-stim power for these events       
        powerData = loadAllPower(yc1Tal,subj,events,params.freqBins,[5100 7000],powParams,eventsToUse,params);
        powerData = permute(powerData,[3 1 2 4]);
        
        % indices of stim and non-stim trials
        stimEvents    = [events(eventsToUse).isStim]==1;
        nonStimEvents = [events(eventsToUse).isStim]==0;
        
        [resTmp,postProb] = computePostStimStats(powerData,yc1Data,chanceData,stimEvents,nonStimEvents,events,eventsToUse,timeToUse,yc1thresh);                
        res(stimLoc).deltaRR         = resTmp.deltaRR;
        res(stimLoc).deltaRR_bin     = resTmp.deltaRR_binary;
        res(stimLoc).yc1Score        = resTmp.yc1Score;
        res(stimLoc).postProbStim    = resTmp.postProbStim;
        res(stimLoc).logOddsStim     = resTmp.logOddsStim;
        res(stimLoc).postProbNonStim = resTmp.postProbNonStim;
        res(stimLoc).logOddsNonStim  = resTmp.logOddsNonStim;
        res(stimLoc).deltaEE         = resTmp.deltaEE;
        res(stimLoc).deltaEE_Prob    = resTmp.deltaEE_Prob;
        res(stimLoc).stimPerf        = resTmp.stimPerf;
        res(stimLoc).nonStimPerf     = resTmp.nonStimPerf;
        res(stimLoc).numNZweights    = resTmp.numNZweights;
        res(stimLoc).AUC             = resTmp.AUC;
                        
        % Do 1000 permutattion of stim non-stim labels
        nIters = 1000;
        stimLabels           = [events(eventsToUse).isStim];
        deltaEE_perm         = NaN(1,nIters);
        deltaEE_Prob_perm    = NaN(1,nIters);
        postProbStim_perm    = NaN(1,nIters);
        postProbNonStim_perm = NaN(1,nIters);
        logOddsStim_perm     = NaN(1,nIters);
        logOddsNonStim_perm  = NaN(1,nIters);
        deltaRR_perm         = NaN(1,nIters);
        deltaRR_bin_perm     = NaN(1,nIters);
        stimPerf_perm        = NaN(1,nIters);
        nonStimPerf_perm     = NaN(1,nIters);
        
        for i = 1:nIters              
            fprintf('Processing %s: stim region: %s, %s (%d of %d). Permuation %d of %d.\n',subj,res(stimLoc).stimAnat,res(stimLoc).stimTagName,stimLoc,length(stimLeads),i,nIters)
            
            if any(strcmp({'first','second'},stimToUse))
                stimLabelsRand  = stimLabels(randperm(length(stimLabels)));
            else
                Ytmp           = reshape(stimLabels,2,[]);
                randOrder      = randperm(length(Ytmp));
                stimLabelsRand = reshape(Ytmp(:,randOrder),[],1);
            end                        
            stimEvents     = stimLabelsRand == 1;
            nonStimEvents  = stimLabelsRand == 0;            
            
            % this is inefficient, don't need to redo glmval each time...
            resTmp = computePostStimStats(powerData,yc1Data,chanceData,stimEvents,nonStimEvents,events,eventsToUse,timeToUse,yc1thresh);              
            deltaEE_perm(i)         = resTmp.deltaEE;
            deltaEE_Prob_perm(i)    = resTmp.deltaEE_Prob;
            deltaRR_perm(i)         = resTmp.deltaRR;
            deltaRR_bin_perm(i)         = resTmp.deltaRR_binary;
            postProbStim_perm(i)    = nanmean(resTmp.postProbStim);
            postProbNonStim_perm(i) = nanmean(resTmp.postProbNonStim);
            logOddsStim_perm(i)     = resTmp.logOddsStim;
            logOddsNonStim_perm(i)  = resTmp.logOddsNonStim;
            stimPerf_perm(i)        = nanmean(resTmp.stimPerf);
            nonStimPerf_perm(i)     = nanmean(resTmp.nonStimPerf);
        end
        res(stimLoc).deltaEE_perm         = deltaEE_perm;
        res(stimLoc).deltaEE_Prob_perm    = deltaEE_Prob_perm;
        res(stimLoc).deltaRR_perm         = deltaRR_perm;
        res(stimLoc).deltaRR_bin_perm     = deltaRR_bin_perm;
        res(stimLoc).postProbStim_perm    = postProbStim_perm;
        res(stimLoc).postProbNonStim_perm = postProbNonStim_perm;
        res(stimLoc).logOddsStim_perm     = logOddsStim_perm;
        res(stimLoc).logOddsNonStim_perm  = logOddsNonStim_perm;         
        res(stimLoc).stimPerf_perm        = stimPerf_perm;
        res(stimLoc).nonStimPerf_perm     = nonStimPerf_perm;
        
        
    end            

    fname = fullfile(saveDir,[subj '_YC2_postStimChange.mat']);
    save(fname,'res');
catch e
    fname = fullfile(saveDir,[subj '_YC2_postStimChange_error.mat']);
    save(fname,'e')
end

function [res,postProb] = computePostStimStats(powerData,yc1Data,chanceData,stimEvents,nonStimEvents,events,eventsToUse,timeToUse,yc1thresh)
% following Youssef's terminology, compute recall relative to
% non-stim baseline. He uses percent recall, here I use performance
% factor change
res = [];
testErrors  = [events(eventsToUse).testError];
res.stimPerf    = 1-testErrors(stimEvents);
res.nonStimPerf = 1-testErrors(nonStimEvents);
res.deltaRR = 100*((mean(res.stimPerf)-mean(res.nonStimPerf))/mean(res.nonStimPerf));

% also make binary
recalled = testErrors<median(testErrors);
% recalled = testErrors<yc1thresh;
res.deltaRR_binary = 100*((mean(recalled(stimEvents))-mean(recalled(nonStimEvents)))/mean(recalled(nonStimEvents)));

% size of feature matrix
nElecs = size(powerData,4);
nFreqs = size(powerData,2);

% Choose most significnat timebin from YC1. Instead of looping
% over all times? If we have a tie, use most sig with highest AUC
if yc1Data.params.doBinary
    field_all = 'auc_all';
    field     = 'AUC';
    link      = 'logit';
else
    field_all = 'r_all';
    field     = 'r';
    link      = 'identity';
end
p = mean(repmat(yc1Data.(field),[size(chanceData.(field_all),1), 1]) > chanceData.(field_all));
AUC = yc1Data.(field);
if strcmpi(timeToUse,'best')    
    maxP = max(p);    
    AUC(p~=maxP) = NaN;
    [~,timeToUse] = max(AUC);
    res.yc1Score = maxP;
else
    timeToUse    = strcmpi(yc1Data.params.timeBinLabels,timeToUse);
    res.yc1Score = p(timeToUse);
end
res.AUC = AUC(timeToUse);

% average weights across all across all YC1 training folds
A                = nanmean(horzcat(yc1Data.res(timeToUse).A{:}),2);
intercept        = nanmean([yc1Data.res(timeToUse).intercept{:}]);
res.numNZweights = nnz(A);

% predict YC1 time bin weights applied to YC2 time bin
B1 = [intercept;A];

% predict with YC1 time bin weights applied to YC2 post period
X = reshape(squeeze(powerData(:,:,1,:)),size(powerData,1),nFreqs*nElecs);
postProb = glmval(B1,X,link);

% convert probability to log-odds, for both stim and non-stim
postProbStim        = postProb(stimEvents);
res.postProbStim    = postProbStim(~isnan(postProbStim));
res.logOddsStim     = nanmean(log(res.postProbStim./(1-res.postProbStim)));
postProbNonStim     = postProb(nonStimEvents);
res.postProbNonStim = postProbNonStim(~isnan(postProbNonStim));
res.logOddsNonStim  = nanmean(log(res.postProbNonStim./(1-res.postProbNonStim)));

% following Youssef's terminology, compute Encoding Efficiency
% change
res.deltaEE      = 100*((res.logOddsStim-res.logOddsNonStim)/abs(res.logOddsNonStim));
res.deltaEE_Prob = 100*((nanmean(res.postProbStim)-nanmean(res.postProbNonStim))/abs(nanmean(res.postProbNonStim)));

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
        fname = fullfile(subjPath,'RAM_YC2_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'_post.mat']);
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




