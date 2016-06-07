function auc = YC2_refactor_classifyStimFromFullModel(subj,modelDir,overwrite)
% function [perf,AUC,subject,params] = YC1_runMulti_subj_ROC(subj,params,saveDir)
%
% Inputs:
%
%       subj - subject string
%     params - params structure
%    savedir - string of path to save directory
%
% Runs lasso regression or classification for one subject using the
% parameters in params. See multiParams for description of parameters.
%
% Saves results to saveDir/<subj>_lasso.mat
%
% Outputs:
%
%     perf: percent classifier accuracy
%      AUC: area under the ROC curve
%  subject: current subject
%   params: params used

% initialize outputs
auc = [];

% do we overwrite results file if it exists?
fname = fullfile(modelDir,[subj '_yc2_class.mat']);
if exist(fname,'file') && ~overwrite
    fprintf('Classifier file already exists for %s.\n',subj)
    return
end

% make sure yc1 model exists
yc1_file = fullfile(modelDir,[subj '_class_pow.mat']);
if ~exist(yc1_file,'file')
    fprintf('No YC1 model for %s. Skipping.\n',subj);
    return
end

% ignore some subejcts with errors
try
    
    % load YC1 model
    yc1Data = load(yc1_file);
    
     % load yc2 events
    if isfield(yc1Data.params,'diffRegress') && yc1Data.params.diffRegress == 1
        events = load(fullfile('/scratch/jfm2/diffRegress/RAM_YC2/',[subj '_events.mat']));
        events = events.events;
    else
        events = get_sub_events('RAM_YC2',subj);    
        % add the test error to the learning trials
        events  = addExtraFields(events);
    end
    
    % filter to just encoding events
    eventsToUse = yc1Data.params.eventFilter(events);    
    
    % if a session is not a full session, exclude it
    sessions  = [events.session];
    evPerSess = grpstats(sessions(eventsToUse),sessions(eventsToUse),'numel');
    uniqSess  = unique(sessions);
    for s = 1:length(uniqSess)
        if evPerSess(s) < 96
            eventsToUse(sessions==uniqSess(s)) = 0;
        end
    end
    
    % filter to just non-stim events
    eventsToUse = eventsToUse & [events.isStim]~=1;
     if sum(eventsToUse) == 0
        fprintf('No enough events for %s.\n',subj)
        return
     end
     
    % new session vector
    sessions  = [events(eventsToUse).session];    
    
    % because we really have only have the number of trials because of the
    % two encoding periods, reduce to just one encoding period
    encPeriods = [events(eventsToUse).withinItemCount]; 
    sessions   = sessions(encPeriods==1);
    
    % response data    
    if isfield(yc1Data.params,'diffRegress') && yc1Data.params.diffRegress == 1
        Y = [events(eventsToUse).residual]';        
        Y = Y > 0;
        Y = Y(encPeriods==1);
    else
        Y = [events(eventsToUse).testError]';        
        Y = Y < median(Y);
        Y = Y(encPeriods==1);
    end    
    
    % load power    
    [powMat,encMat,timeMat,~,tal,powParams] = loadAllPower(yc1Data.tal,subj,events,yc1Data.params.timeBins,eventsToUse,yc1Data.params);    

    % filter to just chosen time and enc
    bestE = yc1Data.res.fullModel.bestE;
    bestT = yc1Data.res.fullModel.bestT;
    featureMask = encMat == bestE & timeMat == bestT;    
    X = powMat(:,featureMask);    
    
    % test on all sessions together, and also each separaetly 
    
    
    % normalize full training set        
    [X,m,s,sessMap] = normPow(X,sessions,[],[],[]);

    % predict on all data
    res     = [];
    res.auc = predict_local(X,Y,yc1Data);

    % get auc null distribution
    res.permAucs = NaN(1,1000);
    for iter = 1:1000 
        randOrder = randperm(length(Y));
        res.permAucs(iter) = predict_local(X,Y(randOrder),yc1Data);
    end        
    fprintf('%s: YC2 AUC = %.3f, p = %.3f.\n',subj,res.auc,mean(res.auc < res.permAucs))
    
    % also predict each YC2 session individually
    uniqSess = unique(sessions);
    res.sess = [];
    for sess = 1:length(uniqSess);
        sessInds = sessions == uniqSess(sess);
        xSess    = X(sessInds,:);
        ySess    = Y(sessInds);
        res.sess(sess).auc = predict_local(xSess,ySess,yc1Data);
        
        % get auc null distribution
        res.sess(sess).permAucs = NaN(1,1000);
        for iter = 1:1000
            randOrder = randperm(length(ySess));
            res.sess(sess).permAucs(iter) = predict_local(xSess,ySess(randOrder),yc1Data);
        end
        fprintf('%s session %d: YC2 AUC = %.3f, p = %.3f.\n',subj,sess,res.sess(sess).auc,mean(res.sess(sess).auc < res.sess(sess).permAucs))
    end    

    if yc1Data.params.saveOutput
        params = yc1Data.params;
        res.bestE = bestE;
        res.bestT = bestT;
        res.bestC = yc1Data.res.fullModel.bestC;
        save(fname,'res','params','tal')
    end
    
catch e
    fname = fullfile(modelDir,[subj '_class_error.mat']);
    save(fname,'e')
end

function auc = predict_local(X,Y,yc1Data)

C = yc1Data.res.fullModel.bestC;
if strcmpi(yc1Data.params.normType,'L1')
    liblin_param = ['-c ' sprintf('%f',C) ' -s 6 -q'];
    trainFun = @train;
    testFun  = @predict;
elseif strcmpi(yc1Data.params.normType,'L2')
    liblin_param = ['-c ' sprintf('%f',C) ' -s 0 -q'];
    trainFun = @train;
    testFun  = @predict;
elseif strcmpi(yc1Data.params.normType,'svm')
    liblin_param = ['-t 2 -g ' sprintf('%f',thisKern) ' -c ' sprintf('%f',C) ' -b 0 -q'];
    trainFun = @svmtrain;
    testFun  = @svmpredict;
end

Y = double(Y);
Y(Y==0) = -1;
[pred, acc, dec] = testFun(double(Y),sparse(X),yc1Data.res.fullModel.model,'-b 1 -q');
[~,~,~,auc]      = perfcurve(Y,dec(:,1),1);

function folds = create_stratKfolds(groups,k)

ind = crossvalind('kfold',groups,k);
folds = true(k,length(groups));
for iFold = 1:k
    folds(iFold,ind==iFold) = false;
end

function [nFolds,folds,trials] = createFolds(session,cvField)

[trials,~,trialInds] = unique([session' cvField'],'rows');
nFolds = size(trials,1);
folds = false(nFolds,size(trialInds,1));
for iFold = 1:nFolds
    folds(iFold,:) = trialInds ~= iFold;
end



% tal = sessData.features.elecs(elecs);
function [powerData,encMat,timeMat,sessInds,tal,powerParams] = loadAllPower(tal,subj,events,timeBins,eventsToUse,params)

% power matrix will initially be # events x # elecs x # times x # freqs
nTimes  = size(timeBins,1);
nEvents = sum(eventsToUse);
[tal]   = filterTalByRegion(tal,params.region);
nElecs  = length(tal);

% load parameters used to create power
powerParams = load(fullfile(params.powerPath,'params_RAM_YC2.mat'));
freqs       = powerParams.params.pow.freqs;
nFreqs      = length(powerParams.params.pow.freqs);
nFreqBins   = size(params.freqBins,1);

% will hold power data
powerDataAllFreqs = NaN(length(events),nElecs,nTimes,nFreqs);
powerData         = NaN(length(events),nElecs,nTimes,nFreqBins);

% 
basePath  = params.powerPath;
subjPath  = fullfile(basePath,subj);
sessions  = [events.session];
uniqSess  = unique(sessions);

% loop over each electrode
for e = 1:nElecs
    fprintf('%s: Loading elec %d of %d.\n',subj,e,nElecs);
    elecNum = tal(e).channel;
    
    % and session
    for s = 1:length(uniqSess)
        sessInds = sessions==uniqSess(s);        
        
        % load power for electrode/session
        fname = fullfile(subjPath,'RAM_YC2_events',num2str(uniqSess(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
        sessPow = load(fname);
        
        % under zscore
        zPow = sessPow.sessOutput.pow;
        nT = size(sessPow.sessOutput.pow,2);
        nE = size(sessPow.sessOutput.pow,3);
        elecPow = zPow .* repmat(sessPow.sessOutput.stdBasePow,[1 nT nE]) + repmat(sessPow.sessOutput.meanBasePow,[1 nT nE]);                
        
        % mean into time bins        
        elecPow = permute(elecPow,[3 2 1]);
        for t = 1:nTimes
            tStart = params.timeBins(t,1);
            tEnd = params.timeBins(t,2);
            tInds = tStart:tEnd;
            powerDataAllFreqs(sessInds,e,t,:) = nanmean(elecPow(:,tInds,:),2);            
        end
    end
end

% bin by frequency
if nFreqBins ~= 0    
    for f = 1:nFreqBins
        fInds = freqs >= params.freqBins(f,1)-.001 & freqs <= params.freqBins(f,2)+.001;
        powerData(:,:,:,f) = nanmean(powerDataAllFreqs(:,:,:,fInds),4);
    end
else
    powerData = powerDataAllFreqs;
end

% make encoding period one and two  a new dimension of matrix
encPeriods = [events.withinItemCount];
powEnc1    = powerData(eventsToUse&encPeriods==1,:,:,:);
powEnc2    = powerData(eventsToUse&encPeriods==2,:,:,:);
powerData  = cat(5,powEnc1,powEnc2);
powerData  = permute(powerData,[1,2,5,3,4]);

% create matrix to keep track of time bins
if nFreqBins == 0;nFreqBins=nFreqs;end
nEncs = size(powerData,3);
timeMat = NaN(1,nElecs,nEncs,nTimes,nFreqBins);
for t = 1:size(timeMat,4)
    timeMat(:,:,:,t,:) = t;
end

% create matrix to keep track of enc periods
encMat = NaN(1,nElecs,nEncs,nTimes,nFreqBins);
for e = 1:size(encMat,3)
    encMat(:,:,e,:,:) = e;
end
 
% reshape into obs x features
powerData = reshape(powerData,size(powerData,1),[]);
timeMat   = reshape(timeMat,[1,nElecs*nEncs*nTimes*nFreqBins]);
encMat    = reshape(encMat,[1,nElecs*nEncs*nTimes*nFreqBins]);

function [powerNorm,m,sigma,sessMap] = normPow(powerData,sessions,m,sigma,sessMap)

[nEvents,nElecs,nEncs,nTimes,nFreqBins] = size(powerData);
powerNorm = NaN(nEvents,nElecs,nEncs,nTimes,nFreqBins);


uniqSess     = unique(sessions);
if isempty(m) || isempty(sigma)
    m       = NaN(length(uniqSess),size(powerData,2));
    sigma   = NaN(length(uniqSess),size(powerData,2));
    sessMap = NaN(length(uniqSess),1);
end

for s = 1:length(uniqSess)    
    sessInds = sessions==uniqSess(s);
    if all(isnan(m(s,:)))
        [powerNorm(sessInds,:),m(s,:),sigma(s,:)] = zscore(powerData(sessInds,:));
        sessMap(s) = uniqSess(s);
    else
        sessMapInd = sessMap == uniqSess(s);
        powerNorm(sessInds,:) = bsxfun(@rdivide, bsxfun(@minus, powerData(sessInds,:), m(sessMapInd,:)), sigma(sessMapInd,:));   
    end
end
    


function events = addExtraFields(events)
% add testError field
% add inner field (1 = inner region, 0 = outer region)

testInd = strcmp({events.type},'NAV_TEST');
recEvents = events(testInd);
[events.testError] = deal(NaN);
[events.recalled] = deal(NaN);
[events.inner] = deal(NaN);
[events.sessionHalf] = deal(0);
[events.withinItemCount] = deal(0);
sessVec = [events.session];
trialVec = [events.blocknum];

uniqSess = unique(sessVec);
sessPerc = NaN(1,length(sessVec));

for sess = 1:length(uniqSess)
    sessInd = sessVec == uniqSess(sess);
    sessPerc(sessInd) = [1:sum(sessInd)]./sum(sessInd);
    
    items = unique(trialVec(sessInd));
    for item = 1:length(items)
        ind = trialVec == items(item) & sessInd;
        
        for c = 1:sum(ind)
            indTmp = find(ind);
            events(indTmp(c)).withinItemCount = c;
        end
    end
    
end
[events(sessPerc>.5).sessionHalf] = deal(1);


for rec = 1:length(recEvents);
    session = recEvents(rec).session;
    trial = recEvents(rec).blocknum;
    err = recEvents(rec).respPerformanceFactor;
    ind = sessVec == session & trialVec == trial;
    [events(ind).testError] = deal(err);
    [events(ind).inner] = deal(abs(recEvents(rec).objLocs(1)) < 568/30 && abs(recEvents(rec).objLocs(2)) < 7);
end














