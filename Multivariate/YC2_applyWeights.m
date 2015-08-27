function [AUC,AUC_enc,perf,perfEnc] = YC2_applyWeights(subj,params,yc1Data,saveDir)
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

AUC     = [];
AUC_enc = [];
perf    = [];
perfEnc = [];

try
    
    % load subject electrode locations and filter to specific regions if
    % desired.
    tal = getBipolarSubjElecs(subj,1,1,1);
    tal = filterTalByRegion(tal,params.region);
       
    % load power parameters    
    powParams = load(fullfile(params.powerPath,'params.mat'));
    
    % Setting time bins for convenience:
    tEnds     = (powParams.params.pow.timeWin:powParams.params.pow.timeStep:powParams.params.eeg.durationMS)+powParams.params.eeg.offsetMS;
    tStarts   = tEnds - powParams.params.pow.timeWin+1;
    powParams.timeBins = [tStarts' tEnds'];
    
    % load events
    events = get_sub_events('RAM_YC2',subj);
    
    % add the test error to the learning trials
    events  = addErrorField(events);
    session = [events.session];
    
    % filter to just NON-STIM learning trials
    eventsToUse = params.eventFilter(events) & [events.isStim]==1; 
    if sum(eventsToUse) < 10
        fprintf('Not enough events for %s.\n',subj)
        return
    end
    
    % get parameters
    freqBins      = params.freqBins;
    timeBins      = params.timeBins; 
    modelEachTime = params.modelEachTime;
    doBinary      = params.doBinary;
    saveOutput    = params.saveOutput;
    doPermute     = params.doPermute;
    
    % load power for all electrodes
    powDir = fullfile(saveDir,'power');
    powFile = fullfile(powDir,[subj '_binnedPower.mat']);
    if ~exist(powDir,'dir')
        mkdir(powDir)
    end
    if params.loadPower
        powerData = load(powFile);
        powerDataEncAvg = powerData.powerDataEncAvg;
        powerData = powerData.powerData;        
    else
        if ~params.useCorrectedPower
            powerData       = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse);
            powerDataEncAvg = loadAllPower(tal,subj,events,freqBins,[1 5000],powParams,eventsToUse);
        else
            powerData       = loadResids_locs(tal,subj,freqBins,timeBins,powParams,eventsToUse);
            powerDataEncAvg = loadResids_locs(tal,subj,freqBins,[1 5000],powParams,eventsToUse);
        end
        powerData       = permute(powerData,[3 1 2 4]);
        powerDataEncAvg = permute(powerDataEncAvg,[3 1 2 4]);
        if params.savePower
            powFile = fullfile(powDir,[subj '_binnedPower.mat']);
            save(powFile,'powerData','powerDataEncAvg','params')
        end
    end
    
    % size of feature matrix
    nElecs = size(powerData,4);
    nTimes = size(powerData,3);
    nFreqs = size(powerData,2);    
    
    % response data
    % SHOULD THIS BE THE MEDIA OF ALL TRIALS OR JUST NON-STIM?
    Y = [events(eventsToUse).testError]';
    if doBinary
        Y  = Y < median(Y);
    end
            
    % permute the responses if desired
    if doPermute
        randOrder = randperm(length(Y)/2);
        Ytmp = reshape(Y,2,[]);
        Y = reshape(Ytmp(:,randOrder),[],1);
    end
    
    objLocs = vertcat(events(eventsToUse).objLocs);
        
    % We can model time points seperately, so # features = # freqs x # elecs,
    % or we can model it all together, so # features = # times x # freqs x #
    % elecs. If we are modeling each time, the weights from the YC1
    % classification at that time bin will be applied to the YC2 data at
    % the same time bin. In addition, the weight will be applied to the
    % average power across the encoding interval.
    res = [];
    if modelEachTime
        perf    = NaN(1,nTimes);
        perfEnc = NaN(1,nTimes);
        for t = 1:nTimes
            
            % average weights across all across all YC1 training folds
            A         = mean(horzcat(yc1Data.res(t).A{:}),2);                  
            intercept = mean([yc1Data.res(t).intercept{:}]);
            
            % reshape power into # trials x # features
            X = reshape(squeeze(powerData(:,:,t,:)),size(powerData,1),nFreqs*nElecs);
            
            % predict YC1 time bin weights applied to YC2 time bin
            B1 = [intercept;A];
            res(t).yPred    = glmval(B1,X,'logit');
            res(t).predBool = res(t).yPred > mean(Y) == Y;
            res(t).perf     = mean(res(t).predBool);
            res(t).A        = A;
            rest(t).intcp   = intercept;
            perf(t)         = res(t).perf;                        
            
            % predict YC1 time bin weights applied to YC2 average encoding
            X_enc = reshape(squeeze(powerDataEncAvg(:,:,1,:)),size(powerDataEncAvg,1),nFreqs*nElecs);
            res(t).yPredEnc    = glmval(B1,X_enc,'logit');
            res(t).predBoolEnc = res(t).yPredEnc > mean(Y) == Y;
            res(t).perfEnc     = mean(res(t).predBoolEnc);                        
            perfEnc(t)         = res(t).perfEnc;
                      
            % calculate area under ROC curve
            if doBinary                                          
                [res(t).xAUC,res(t).yAUC,~,res(t).AUC,res(t).optPoint] = perfcurve(Y,res(t).yPred,true);
                AUC(t) = res(t).AUC;
                
                [res(t).xAUC_enc,res(t).yAUC_enc,~,res(t).AUC_enc,res(t).optPoint_enc] = perfcurve(Y,res(t).yPredEnc,true);
                AUC_enc(t) = res(t).AUC_enc;
            end                        
        end        

        % if using all time points in one model
    else
        
        % average weights across all across all YC1 training folds
        A         = mean(horzcat(yc1Data.res.A{:}),2);
        intercept = mean([yc1Data.res.intercept{:}]);
        
        % reshape into # trials x # features
        X = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs);
                
        % predict YC1 time bin weights applied to YC2 time bin
        B1 = [intercept;A];
        res.yPred    = glmval(B1,X,'logit');
        res.predBool = res.yPred > mean(Y) == Y;
        res.perf     = mean(res.predBool);
        res.A        = A;
        res.intcp    = intercept;
        perf         = res.perf;
                       
        % calculate area under ROC curve
        if doBinary
            [res.xAUC,res.yAUC,~,res.AUC,res.optPoint] = perfcurve(Y,res.yPred,true);
            AUC = res.AUC;
        end        
    end
        
    subject       = subj;    
    if saveOutput
        fname = fullfile(saveDir,[subj '_YC2_lasso.mat']);
        save(fname,'res','Y','objLocs','params','perf','tal','AUC','AUC_enc','perfEnc');
    end
catch e
    fname = fullfile(saveDir,[subj '_YC2_lasso_error.mat']);
    save(fname,'e')
end

function powerData = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse)

nFreqs = size(freqBins,1);
nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
nElecs = length(tal);
powerData = NaN(nFreqs,nTimes,nEvents,nElecs);

for e = 1:nElecs
    elecNum = tal(e).channel;
      
    basePath  = '/data10/scratch/jfm2/RAM/biomarker/power/';
    subjPath  = fullfile(basePath,subj);
    sessions = unique([events.session]);
    subjPow  = [];
    for s = 1:length(sessions)
       fname = fullfile(subjPath,'RAM_YC2_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
       sessPow = load(fname);
       subjPow = cat(3,subjPow,sessPow.sessOutput.pow);
    end
    
    if length(eventsToUse) ~= size(subjPow,3)
        fprintf('Number of events does not match size of power matrix for %s!.\n',subj)
        return
    end    
    subjPow = subjPow(:,:,eventsToUse);
    
    % average frequencies
    tmpPower = NaN(nFreqs,size(subjPow,2),size(subjPow,3));
    for f = 1:nFreqs
        fInds = powParams.params.pow.freqs >= freqBins(f,1) & powParams.params.pow.freqs < freqBins(f,2);
        tmpPower(f,:,:) = nanmean(subjPow(fInds,:,:),1);
    end
    subjPow = tmpPower;
    
    % average times
    tmpPower = NaN(nFreqs,nTimes,size(subjPow,3));
    for t = 1:nTimes
        tInds = powParams.timeBins(:,1) >= timeBins(t,1) & powParams.timeBins(:,2) < timeBins(t,2);
        tmpPower(:,t,:) = nanmean(subjPow(:,tInds,:),2);
    end
    powerData(:,:,:,e) = tmpPower;
end

function [powerData] = loadResids_locs(tal,subj,freqBins,timeBins,powParams,eventsToUse)

nFreqs = size(freqBins,1);
nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
nElecs = length(tal);
powerData = NaN(nFreqs,nTimes,nEvents,nElecs);

basePath  = '/data10/scratch/jfm2/YC1/multi/power/regress/';
subjPath  = fullfile(basePath,subj);

for e = 1:nElecs
    elecNum   = tal(e).channel;
    fname     = sprintf('%s_elec_%d-%d_residuals.mat',subj,elecNum(1),elecNum(2));
    
    if ~exist(fullfile(subjPath,fname),'file')
        error('Residuals file %s not found.\n',fname)
    else
        elecData = load(fullfile(subjPath,fname));
        resids   = elecData.resid;
        resids   = permute(resids,[3 2 1]);
        if size(resids,3) ~= sum(eventsToUse)
            keyboard
        end
        
        % average frequencies
        tmpPower = NaN(nFreqs,size(resids,2),size(resids,3));
        for f = 1:nFreqs
            fInds = powParams.params.pow.freqs >= freqBins(f,1) & powParams.params.pow.freqs < freqBins(f,2);
            tmpPower(f,:,:) = nanmean(resids(fInds,:,:),1);
        end
        resids = tmpPower;
        
        % average times
        tmpPower = NaN(nFreqs,nTimes,size(resids,3));
        for t = 1:nTimes
            tInds = powParams.timeBins(:,1) >= timeBins(t,1) & powParams.timeBins(:,2) < timeBins(t,2);
            tmpPower(:,t,:) = nanmean(resids(:,tInds,:),2);
        end
        powerData(:,:,:,e) = tmpPower;
        
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




