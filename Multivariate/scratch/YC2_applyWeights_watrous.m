function [AUC,AUC_enc,perf,perfEnc,AUC_best,AUC_bestEnc] = YC2_applyWeights(subj,params,yc1Data,yc1ChanceData,saveDir)
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

AUC         = [];
AUC_enc     = [];
perf        = [];
perfEnc     = [];
AUC_best    = [];
AUC_bestEnc = [];
try
    
    % load subject electrode locations and filter to specific regions if
    % desired.
    tal = getBipolarSubjElecs(subj,1,1,1);
    tal = filterTalByRegion(tal,params.region);
    if isempty(tal)
        fprintf('No %s electrode for %s.\n',params.region,subj)
        return
    end
    
    % load power parameters    
    % WE ARE PROBABLY GOING TO NEED SEPERATE YC1 AND YC2 params
    powParams = load(fullfile(params.powerPath,'params.mat'));
    
    % Setting time bins for convenience:
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
        powerData       = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params);
        powerDataEncAvg = loadAllPower(tal,subj,events,freqBins,[1 5000],powParams,eventsToUse,params);
        powerData       = permute(powerData,[3 1 2 4]);
        powerDataEncAvg = permute(powerDataEncAvg,[3 1 2 4]);
        
        % load andrew's data
        net = load(fullfile('/scratch/jfm2/network',[subj '_Sess0_Results.mat']));
        wstPPC_size = size(net.results.weighted_stPPC);
        wstPPC = NaN(wstPPC_size(1),4,6,wstPPC_size(4));
        tbins = {[1:20],[21:40],[41:60],[61:80],[81:100],[1:100]};
        fbins = {[1:3],[4:8],[11:16],[18:22]};
        
        for f = 1:4
            tmp = nanmean(net.results.weighted_stPPC(:,fbins{f},:,:),2);
            for t = 1:6
                wstPPC(:,f,t,:) = nanmean(tmp(:,:,tbins{t},:),3);
            end
        end
        
        wstPPC = permute(wstPPC,[4 2 3 1]);
        powerData = cat(4,powerData,wstPPC);
        powerDataEncAvg = cat(4,powerDataEncAvg,wstPPC(:,:,6,:));
        
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
    % OR the median of the YC1 data?
    Y = [events(eventsToUse).testError]';

    if doBinary
        Y  = Y < yc1thresh;
    end
            
    % permute the responses if desired
    if doPermute
        randOrder = randperm(length(Y)/2);
        Ytmp = reshape(Y,2,[]);
        Y = reshape(Ytmp(:,randOrder),[],1);
    end
    
    objLocs = vertcat(events(eventsToUse).objLocs);
    
    % What was the best YC1 time point?
    p = mean(repmat(yc1Data.AUC,[size(yc1ChanceData.auc_all,1), 1]) > yc1ChanceData.auc_all);
    maxP = max(p);
    yc1AUC = yc1Data.AUC;
    yc1AUC(p~=maxP) = NaN;
    [~,bestTime] = max(yc1AUC);
    
    
    % We can model time points seperately, so # features = # freqs x # elecs,
    % or we can model it all together, so # features = # times x # freqs x #
    % elecs. If we are modeling each time, the weights from the YC1
    % classification at that time bin will be applied to the YC2 data at
    % the same time bin. In addition, the weight will be applied to the
    % average power across the encoding interval.
    res     = [];
    resBest = [];
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
                [~,~,~,res(t).AUC] = perfcurve(Y,res(t).yPred,true);
                AUC(t) = res(t).AUC;
                
                [~,~,~,res(t).AUC_enc] = perfcurve(Y,res(t).yPredEnc,true);
                AUC_enc(t) = res(t).AUC_enc;
            end
            
            % average weights across all across all YC1 training folds for the
            % best timepoint
            A         = mean(horzcat(yc1Data.res(bestTime).A{:}),2);
            intercept = mean([yc1Data.res(bestTime).intercept{:}]);
                        
            % predict YC1 best bin weights applied to YC2 time bin
            B1 = [intercept;A];
            resBest(t).yPred    = glmval(B1,X,'logit');
            resBest(t).predBool = resBest(t).yPred > mean(Y) == Y;
            resBest(t).perf     = mean(resBest(t).predBool);
            resBest(t).A        = A;
            resBest(t).intcp    = intercept;
            
            % predict YC1 time bin weights applied to YC2 average encoding            
            resBest(t).yPredEnc    = glmval(B1,X_enc,'logit');
            resBest(t).predBoolEnc = resBest(t).yPredEnc > mean(Y) == Y;
            resBest(t).perfEnc     = mean(resBest(t).predBoolEnc);
            
            % calculate area under ROC curve
            if doBinary
                [~,~,~,resBest(t).AUC] = perfcurve(Y,resBest(t).yPred,true);
                AUC_best(t) = resBest(t).AUC;
                
                [~,~,~,resBest(t).AUC_enc] = perfcurve(Y,resBest(t).yPredEnc,true);
                AUC_bestEnc(t) = resBest(t).AUC_enc;
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
            [~,~,~,res.AUC] = perfcurve(Y,res.yPred,true);
            AUC = res.AUC;
        end        
    end
        
    subject       = subj;    
    if saveOutput
        fname = fullfile(saveDir,[subj '_YC2_lasso.mat']);
        save(fname,'res','Y','objLocs','params','perf','tal','AUC','AUC_enc','perfEnc','resBest','AUC_best','AUC_bestEnc');
    end
catch e
    fname = fullfile(saveDir,[subj '_YC2_lasso_error.mat']);
    save(fname,'e')
end

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
      
    basePath  = '/data10/scratch/jfm2/RAM/biomarker/power/';
    subjPath  = fullfile(basePath,subj);
    sessions = unique([events.session]);
    subjPow  = [];
    for s = 1:length(sessions)
       fname = fullfile(subjPath,'RAM_YC2_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
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
    tmpPower = NaN(nFreqs,nTimes,size(subjPow,3));
    for t = 1:nTimes
        tInds = powParams.timeBins(:,1) >= timeBins(t,1) & powParams.timeBins(:,2) < timeBins(t,2);
        tmpPower(:,t,:) = nanmean(subjPow(:,tInds,:),2);
    end
    powerData(:,:,:,e) = tmpPower;
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




