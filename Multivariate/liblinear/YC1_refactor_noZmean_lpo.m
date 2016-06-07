function auc = YC1_refactor_noZmean_lpo(subj,params,saveDir)
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
fname = fullfile(saveDir,[subj '_class_pow.mat']);
if params.usePhase==1
    fname = fullfile(saveDir,[subj '_class_phase.mat']);
elseif params.usePhase==2
    fname = fullfile(saveDir,[subj '_class_powphase.mat']);
end
if exist(fname,'file') && params.saveOutput && ~params.overwrite
    fprintf('Classifier file already exists for %s.\n',subj)
    return
end

% Using Andrew's power and phase values.
featureDir = fullfile('/home2/andrew.watrous/Results/YC1_Features',subj,'Features_4bands_all_elecs_11_16_2015');
if ~exist(featureDir,'dir') && params.useWatrous
    fprintf('no data for %s.\n',subj)
    return
end

% ignore some subejcts with errors
try
    
    % get parameters    
    saveOutput    = params.saveOutput;
    doPermute     = params.doPermute;
    normType      = params.normType;
    useKfold      = params.useKfold;
    cvField       = params.cvField;
    nestedCvField = params.nestedCvField;
    prctileThresh = params.auc_prctileThresh;
    Cs            = params.Cs;
    Gs            = params.Gs;
    
     % load events
    if isfield(params,'diffRegress') && params.diffRegress == 1
        events = load(fullfile('/scratch/jfm2/diffRegress/RAM_YC1/',[subj '_events.mat']));
        events = events.events;
    else
        events = get_sub_events('RAM_YC1',subj);    
        % add the test error to the learning trials
        events  = addExtraFields(events);
    end
    
    % filter to just encoding events
    eventsToUse = params.eventFilter(events);
    
    % vectors for sessions,and selector
    sessions   = [events(eventsToUse).session]; 
    sel        = [events(eventsToUse).(cvField)];
    
    % because we really have only have the number of trials because of the
    % two encoding periods, reduce to just one encoding period
    encPeriods = [events(eventsToUse).withinItemCount]; 
    sessions   = sessions(encPeriods==1);
    sel        = sel(encPeriods==1);
    
    % response data    
    if isfield(params,'diffRegress') && params.diffRegress == 1
        Y = [events(eventsToUse).residual]';        
        Y = Y > 0;
        Y = Y(encPeriods==1);
    else
        Y = [events(eventsToUse).testError]';        
        Y = Y < median(Y);
        Y = Y(encPeriods==1);
    end    
    
    % load power
    tal = getBipolarSubjElecs(subj,1,1,params.excludeEpiElecs);
    [powMat,encMat,timeMat,~,tal,powParams] = loadAllPower(tal,subj,events,params.timeBins,eventsToUse,params);    
%     nElecs = size(powMat,2);
%     nEncs  = size(powMat,3);
%     nTimes = size(powMat,4);
%     nFreqs = size(powMat,5);
%     

    
    % reshape time and enc matrices into # obs x features. We'll reshape
    % the power values later, since we are to normalize and average the
    % frequencies seperately for each train/test split
%     X = reshape(powMat,size(powMat,1),[]);
%     T = reshape(timeMat,[1,nElecs*nEncs*nTimes*nFreqs]);
%     E = reshape(timeMat,[1,nElecs*nEncs*nTimes*nFreqs]);
    
    
    
%     X = load('/scratch/jfm2/python/YC1/RAM_YC1_R1001P/pow.mat')
%     X = X.pow_mat;
%     T = load('/scratch/jfm2/python/YC1/RAM_YC1_R1001P/time.mat')
%     T = T.time_mat;
    
   

    % create cross val folds
    [nFolds,folds,trials] = createFolds(sessions,sel);
    
    yProbs = [];
    yTests = [];
    for thisFold = 1:nFolds
        fprintf('%s: fold %d of %d.\n',subj,thisFold,nFolds)
        [yProbFold,yTestFold] = run_fold_with_nested_folds(powMat,encMat,timeMat,Y,sessions,folds,thisFold,sel,params,powParams);        
        yProbs = [yProbs,yProbFold];
        yTests = [yTests,yTestFold];
        if ismember([-1,1],yTests)
            [~,~,~,aucTmp] = perfcurve(yTests,yProbs,1);
            fprintf('%s: running AUC = %.3f\n',subj,aucTmp)
        end
       
    end
    
    auc = aucTmp;
    save(fname,'yProbs','yTests','auc','params')
    

catch e
    fname = fullfile(saveDir,[subj '_class_error.mat']);
    save(fname,'e')
end

function [yProb,yTest] = run_fold_with_nested_folds(X,E,T,Y,sessions,folds,thisFold,sel,params,powParams)

% convert Y into double for liblinear
Y = double(Y);
Y(Y==0) = -1;

% unique times and enc periods
Ts = 1:size(params.timeBins,1);
Es = params.encPeriods;
Cs = params.Cs;

% get training data
trainMask = folds(thisFold,:);
[~,subFolds] = createFolds(sessions(trainMask),sel(trainMask));
xTrain    = X(trainMask,:);
yTrain    = Y(trainMask);
sessTrain = sessions(trainMask);

% run nested fold to find optimal enc, time bin, and penalty
[aucs,toRemove] = run_nested_fold(xTrain,E,T,yTrain,sessTrain,subFolds,params,powParams);
ind=find(aucs==max(aucs(:)),1,'first');
[eBestInd,tBestInd,cBestInd]=ind2sub(size(aucs),ind);
thisE = Es(eBestInd);
thisT = Ts(tBestInd);
thisC = Cs(cBestInd);
fprintf('Max AUC = %.3f, best enc = %d, best time = %d, best C = %.4f\n',max(aucs(:)),thisE,thisT,thisC)

% balance the classes
xTrain(toRemove,:) = [];
yTrain(toRemove) = [];
sessTrain(toRemove) = [];

% normalize full training set        
[xTrain,m,s,sessMap] = normPow(xTrain,sessTrain,[],[],[]);
% [xTrain,T,E,m,s,sessMap] = normPow_andMeanFreqs(xTrain,sessTrain,params,powParams.params.pow.freqs,[],[],[]);

% mask encoding period and time bin
featureMask = E == thisE & T == thisT;
xTrain = xTrain(:,featureMask);

% set liblinear inputs
if strcmpi(params.normType,'L1')
    liblin_param = ['-c ' sprintf('%f',thisC) ' -s 6 -q'];
    trainFun = @train;
    testFun  = @predict;
elseif strcmpi(params.normType,'L2')
    liblin_param = ['-c ' sprintf('%f',thisC) ' -s 0 -q'];
    trainFun = @train;
    testFun  = @predict;
elseif strcmpi(params.normType,'svm')
    liblin_param = ['-t 2 -g ' sprintf('%f',thisKern) ' -c ' sprintf('%f',thisC) ' -b 0 -q'];
    trainFun = @svmtrain;
    testFun  = @svmpredict;
end

% create model
% weight by percentage of positive and negative classes
% pos = mean(yTrain==1);
% neg = mean(yTrain==-1);
% mean_tmp = mean([pos neg]);
% pos = sprintf('%f',pos/mean_tmp);
% neg = sprintf('%f',neg/mean_tmp);
% param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
param = liblin_param;
model = trainFun(double(yTrain),sparse(xTrain),param);

% test data
xTest    = X(~trainMask,:);
yTest    = Y(~trainMask);
sessTest = sessions(~trainMask);

% normalize test data and mean into freqs of interest
xTest = normPow(xTest,sessTest,m,s,sessMap);
% xTest = normPow_andMeanFreqs(xTest,sessTest,params,powParams.params.pow.freqs,m,s,sessMap);

% test on held out
xTest = xTest(:,featureMask);
[pred, acc, dec] = testFun(double(yTest),sparse(xTest),model,'-b 1 -q');

% probability of positive class
yProb = dec(:,model.Label==1);

function [aucs,toRemove] = run_nested_fold(X,E,T,Y,sessions,subFolds,params,powParams)

% convert Y into double for liblinear
Y = double(Y);
Y(Y==0) = -1;

% balance the classes
pos = find(Y==1);
neg = find(Y==-1);
numToRemove = length(pos)-length(neg);
toRemove = [];
if numToRemove > 0
    toRemove=randsample(pos,abs(numToRemove));
elseif numToRemove < 0
    toRemove=randsample(neg,abs(numToRemove));
end

Y(toRemove) = [];
X(toRemove,:) = [];
% subFolds(toRemove,:) = [];
% subFolds(:,toRemove) = [];

% find all pairs
pos = find(Y==1);
neg = find(Y==-1);
[A,B] = meshgrid(pos,neg);
pairs = [A(:) B(:)];

% unique times and enc periods
Ts = 1:size(params.timeBins,1);
Es = params.encPeriods;
Cs = params.Cs;

% will hold pPank for each combination of parametres
pRank = NaN(length(Es),length(Ts),length(Cs),length(pairs));

% max loop number for parfor
maxInd = length(Cs)*length(Ts)*length(Es);

parfor thisP = 1:size(pairs,1)
    
    % training data with pair held out
    obsMask = true(1,length(Y));
    obsMask(pairs(thisP,:)) = false;
    
    % training data for fold
    xTrain    = X(obsMask,:);
    yTrain    = Y(obsMask);
    sessTrain = sessions(obsMask);
    
    % normalize training data and mean into frequencies of interest
    [xTrain,m,s,sessMap] = normPow(xTrain,sessTrain,[],[],[]);
    
    % test data
    xTest    = X(~obsMask,:);
    yTest    = Y(~obsMask);
    sessTest = sessions(~obsMask);
    
    % normalize test data
    xTest = normPow(xTest,sessTest,m,s,sessMap);
    
    pRank_fold = NaN(length(Es),length(Ts),length(Cs));
    for ind = 1:maxInd
        
        % get current, E, T, C
        [thisE,thisT,thisC] = ind2sub([length(Es),length(Ts),length(Cs)],ind);
        thisE = Es(thisE);
        thisT = Ts(thisT);
        thisC = Cs(thisC);
        
        % set liblinear inputs
        if strcmpi(params.normType,'L1')
            liblin_param = ['-c ' sprintf('%f',thisC) ' -s 6 -q'];
            trainFun = @train;
            testFun  = @predict;
        elseif strcmpi(params.normType,'L2')
            liblin_param = ['-c ' sprintf('%f',thisC) ' -s 0 -q'];
            trainFun = @train;
            testFun  = @predict;
        elseif strcmpi(params.normType,'svm')
            liblin_param = ['-t 2 -g ' sprintf('%f',thisKern) ' -c ' sprintf('%f',thisC) ' -b 0 -q'];
            trainFun = @svmtrain;
            testFun  = @svmpredict;
        end                
        
        % mask encoding period and time bin
        featureMask = E == thisE & T == thisT;
        xTrain_fold = xTrain(:,featureMask);
        
        % create model
        param = liblin_param;
        model = trainFun(double(yTrain),sparse(xTrain_fold),param);
        
        % test on held out
        xTest_fold = xTest(:,featureMask);
        [pred, acc, dec] = testFun(double(yTest),sparse(xTest_fold),model,'-b 1 -q');        
        pRank_fold(ind) = dec(yTest==1,1) > dec(yTest==-1,1);                
    end
    pRank(:,:,:,thisP) = pRank_fold;
%     ind = sub2ind(size(pRank),[e,t,c,thisP])
%     
%     dec_values = [];
%     labels = [];
%     for thisFold = 1:size(subFolds,1)
%         
%         
%         
%         dec_values = [dec_values; dec(:,model.Label==1)];
%         labels = [labels;yTest];        
%     end  
%     [~,~,~,aucs(ind)] = perfcurve(labels,dec_values,1);
end
aucs = mean(pRank,4);


function [nFolds,folds,ind] = createKfolds(n,percentCV)
if percentCV > .5
    percentCV = .5;
end
ind = crossvalind('kfold',n,round(n/(n*percentCV)));
nFolds = length(unique(ind));
folds = true(nFolds,n);
for iFold = 1:nFolds
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
powerParams = load(fullfile(params.powerPath,'params_RAM_YC1.mat'));
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
        fname = fullfile(subjPath,'RAM_YC1_events',num2str(uniqSess(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
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














