function [perf,AUC,subject,params] = YC1_watrous_phase_wSubSamp(subj,params,saveDir)
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
perf    = [];
subject = [];
AUC     = [];

% do we overwrite results file if it exists?
fname = fullfile(saveDir,[subj '_lasso_pow.mat']);
if params.usePhase==1
    fname = fullfile(saveDir,[subj '_lasso_phase.mat']);
elseif params.usePhase==2
    fname = fullfile(saveDir,[subj '_lasso_powphase.mat']);
end
if exist(fname,'file') && params.saveOutput && ~params.overwrite
    fprintf('Lasso file already exists for %s.\n',subj)
    return
end

% Using Andrew's power and phase values.
featureDir = fullfile('/home2/andrew.watrous/Results/YC1_Features',subj,'Features_4bands_all_elecs_11_16_2015');
if ~exist(featureDir,'dir')
    fprintf('no data for %s.\n',subj)
    return
end

% ignore some subejcts with errors
try
    
    % get parameters
    doBinary      = params.doBinary;
    saveOutput    = params.saveOutput;
    doPermute     = params.doPermute;
    normType      = params.normType;
    useKfold      = params.useKfold;
    cvField       = params.cvField;
    nestedCvField = params.nestedCvField;
    prctileThresh = params.auc_prctileThresh;
    Cs            = params.Cs;
    
    % if we are given one penality parameter and one time bin, and we are
    % not cross validating the encoding period, then we don't need an inner
    % CV loop
    tBest = [];
    encBest = [];
    lambda = [];
    doCalcPenalty = 1;
    if length(Cs) == 1 && size(params.timeBins,1) == 1 && ~strcmp(params.encPeriod,'both')
        doCalcPenalty = 0;
        lambda = Cs;
        tBest = 1;
        if strcmp(params.encPeriod,'first')
            encBest = 1;
        elseif strcmp(params.encPeriod,'second')            
            encBest = 2;
        elseif strcmp(params.encPeriod,'combined')            
            encBest = 3;
        end
    end
    
    % load events
    events = get_sub_events('RAM_YC1',subj);
    
    % add the test error to the learning trials
    events  = addExtraFields(events);
    session = [events.session];
    
    % SOME CHECKS TO ENSURE ENOUGH EVENTS AND SESSIONS
    %---------------------------------------------------------------------%
    eventsToUse = params.eventFilter(events);
    if sum(eventsToUse) < 10
        fprintf('Not enough events for %s.\n',subj)
        return
    end
    encPeriod = [events(eventsToUse).withinItemCount];
    
    if strcmpi(nestedCvField,'session')
        if (length(unique(session)) <= 2 && doCalcPenalty) || (~doCalcPenalty && length(unique(session)) == 1)
            fprintf('Not enough sessions for session nested crossval for %s.\n',subj)
            return
        end
    end
    if strcmpi(cvField,'session')
        if length(unique(session)) < 2
            fprintf('Not enough sessions for session outer crossval for %s.\n',subj)
            return
        end
    end
    %---------------------------------------------------------------------%
    
    % load power or phase
    if params.useWatrous == 1
        [powerData,sessInds,~,tal] = loadWatrousFeatures(subj,params);
    else
        tal = getBipolarSubjElecs(subj,1,1,params.excludeEpiElecs);
        [powerData,sessInds] = loadAllPower(tal,subj,events,params.timeBins,eventsToUse,params);
        powerData = permute(powerData,[3 1 2 4]);
    end
    
    if isempty(tal)
        fprintf('no mtl elecs for %s.\n',subj)
        return
    end
    
    % make sure the power values correspond to the events I expect
    if sum(eventsToUse) ~= size(powerData,1)
        fprintf('EVENTS MISMATCH FOR %s.\n',subj)
        return
    end
    
    % matrix to keep track of encoding trial 1 or 2 feature
    firstEnc  = powerData(encPeriod==1,:,:,:);
    secondEnc = powerData(encPeriod==2,:,:,:);
    
    if strcmp(params.encPeriod,'first')
        powerData = firstEnc;
        encLabel  = ones(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5));
    elseif strcmp(params.encPeriod,'second')
        powerData = secondEnc;
        encLabel  = ones(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5))+1;
    elseif strcmp(params.encPeriod,'both')
        powerData = cat(5,firstEnc,secondEnc);
        encLabel = NaN(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5));
        for t = 1:2
            encLabel(:,:,:,:,t) = t;
        end
    elseif strcmp(params.encPeriod,'combined')
        powerData = cat(5,firstEnc,secondEnc);
        encLabel = ones(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5))+2;
    end
    
    % create a matrix to keep track of the time bin features
    timeLabel = NaN(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5));
    for t = 1:size(timeLabel,3)
        timeLabel(:,:,t,:,:) = t;
    end
    
    % size of feature matrix
    nElecs = size(powerData,4);
    nTimes = size(powerData,3);
    nFreqs = size(powerData,2);
    nItemObs = size(powerData,5);
    
    % response data
    sessInds = sessInds(encPeriod==1);
    Y = [events(eventsToUse).testError]';
    Y = Y(encPeriod==1);
    if doBinary
        Y  = Y < median(Y);
    end
    
    % determine the cross validation folds.
    session = session(eventsToUse);
    cvGroup = [events(eventsToUse).(cvField)];
    nestedCvGroup = [events(eventsToUse).(nestedCvField)];
    nestedCvGroup = nestedCvGroup(encPeriod==1);
    
    if ~useKfold
        [nFolds,folds,trials] = createFolds(session(encPeriod==1),cvGroup(encPeriod==1));
    else
        [nFolds,folds,trials] = createKfolds(sum(encPeriod==1),params.percentCV);
    end
    
    % permute the responses if desired
    if doPermute
        randOrder = randperm(size(trials,1));
        if ~isfield(params,'encPeriod') || strcmpi(params.encPeriod,'both')
            randOrder = [randOrder;randOrder];
            randOrder = randOrder(:);
            Y = Y(randOrder);
        else
            Y = Y(randOrder);
        end
    end
    
    
    % reshape into obs x features
    res = [];
    X         = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs*nItemObs);
    T         = reshape(squeeze(timeLabel),size(timeLabel,1),nFreqs*nTimes*nElecs*nItemObs);
    trialType = reshape(squeeze(encLabel),size(encLabel,1),nFreqs*nTimes*nElecs*nItemObs);
    
 
    % run for each fold
    [res.yProb,res.yPred,res.yTest,res.A,res.intercept,res.err,...
        res.lambda,res.tBest,res.encBest,res.Cs,res.Ts,res.aucs] = deal(cell(nFolds,1));
    for iFold = 1:nFolds
        fprintf('Subject %s: Fold %d of %d.\n',subj,iFold,nFolds)
        [res.yProb{iFold},...
            res.yPred{iFold},...
            res.yTest{iFold},...
            res.A{iFold},...
            res.err{iFold},...
            res.lambda{iFold},...
            res.tBest{iFold},...
            res.encBest{iFold},...
            res.Cs{iFold},...
            res.Ts{iFold},...
            res.aucs{iFold}] = doRegFun(X,Y,T,folds,iFold,lambda,normType,sessInds,trialType,prctileThresh,nestedCvGroup,params.percentCV,useKfold,Cs,tBest,encBest,doCalcPenalty);
        fprintf('Subject %s: Fold %d of %d, %.3f running percent correct.\n',subj,iFold,nFolds,mean(vertcat(res.err{:})))
        
    end
    
    yTest = vertcat(res.yTest{:});
    yProb = vertcat(res.yProb{:});
    [~,~,~,res.AUC] = perfcurve(yTest,yProb,1);
    AUC = res.AUC;
    perf = mean(vertcat(res.err{:}));
    
    % Finally, create a model based on the full dataset for the patient. If
    % modelEachTime, do it for the optimal time
    %     bestTime = find(AUC == max(AUC),1,'first');
    %     X = reshape(squeeze(powerData(:,:,bestTime,:,:)),size(powerData,1),nFreqs*nElecs*nItemObs);
    %     resFull = doRegFullData(X,Y,T,folds,normType,sessInds);
    
    subject       = subj;
    params.lambda = lambda;
    if saveOutput
        objLocs = vertcat(events(eventsToUse).objLocs);
        save(fname,'res','Y','objLocs','params','perf','tal','AUC');
    end
catch e
    fname = fullfile(saveDir,[subj '_lasso_error.mat']);
    save(fname,'e')
end


function [cBest,tBest,encBest,Cs,Ts,auc_pen] = calcPenalty(X,Y,T,folds,sessInds,trialType,normType,prctileThresh,Cs)
% Calculate the optimal penalty parameter.
%
% Returns: penalty (the optimal C)
%          C (vector C values tried)
%          auc_pen (the auc for each penalty paramter)

Y = double(Y);
Y(Y==0) = -1;
if isempty(Cs)
    if strcmpi(normType,'L1')
        Cs = logspace(log10(1e-2),log10(1e4),22);
    elseif strcmpi(normType,'L2')
        Cs = logspace(log10(1e-6),log10(1e4),22);
    end
end
Ts      = 1:length(unique(T));
nCV     = size(folds,1);
types   = unique(trialType);
auc_pen = NaN(length(Cs),length(Ts),length(types));
% pcorr   = NaN(length(Cs),length(Ts),length(types));

% dec_values = NaN(length(Cs),length(Ts),length(types),nCV);
% labels     = NaN(length(Cs),length(Ts),length(types),nCV);
% preds      = NaN(length(Cs),length(Ts),length(types),nCV);

maxInd = length(Cs)*length(Ts)*(length(types));
% maxInd = length(Cs)*length(Ts)*(length(types)*nCV);
 
parfor ind = 1:maxInd
    
    [thisC,thisT,thisEnc] = ind2sub([length(Cs),length(Ts),length(types)],ind);
%     [thisC,thisT,thisEnc,cv] = ind2sub([length(Cs),length(Ts),length(types),nCV],ind);
    thisPen = Cs(thisC);
    
    % set parameters for either L1 or L2 with the current c
    if strcmpi(normType,'L1')
        liblin_param = ['-c ' sprintf('%f',thisPen) ' -s 6 -q'];
    elseif strcmpi(normType,'L2')
        liblin_param = ['-c ' sprintf('%f',thisPen) ' -s 0 -q'];
    end
    
    tInds = T==Ts(thisT);
    if thisEnc > length(types)
        encInds = tInds;
    else
        encInds = trialType==types(thisEnc);
    end
    
    % hold test output
    dec_values = [];%NaN(nCV,1);
    labels = [];%NaN(nCV,1);
    preds = [];%NaN(nCV,1);
    
%     models = [];
    for cv = 1:nCV
        
        inds = folds(cv,:);
        
        % train data for this cv
        xTrain = X(inds,tInds&encInds);
        %         xTrain = [ones(sum(inds),1) xTrain];
        yTrain = Y(inds);
        [xTrain,m,s] = standardize(xTrain,sessInds(inds));
        
        % test data for this cv
        xTest =  X(~inds,tInds&encInds);
        %         xTest = [ones(sum(~inds),1) X(~inds,:)];
        yTest = Y(~inds);
        xTest = standardize_test(xTest,sessInds(~inds),m,s);
        
        
        % weight by percentage of positive and negative classes. I don't
        % totally get this but I am following Jim's code
        pos = mean(yTrain==1);
        neg = mean(yTrain==-1);
        mean_tmp = mean([pos neg]);
        pos = sprintf('%f',pos/mean_tmp);
        neg = sprintf('%f',neg/mean_tmp);
        param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
        
        % train on this data with this cv
        model = train(double(yTrain),sparse(xTrain),param);
%         models = [models;model.w];
        
        % predict with model
        [pred, acc, dec] = predict(double(yTest),sparse(xTest),model,'-b 1 -q');
%         dec_values(ind) = dec(:,1);
%         labels(ind) = yTest;
%         preds(ind) = pred;
        dec_values = [dec_values; dec(:,1)];
        labels = [labels;yTest];
%         preds = [preds;pred];
%         keyboard
    end
    
    [~,~,~,auc_pen(ind)] = perfcurve(labels,dec_values,1);
%     pcorr(ind) = mean(preds == labels);
    
end

% for c = 1:size(dec_values,1)
%     for t = 1:size(dec_values,2)
%         for e = 1:size(dec_values,3);
%            [~,~,~,auc_pen(c,t,e)] = perfcurve(squeeze(labels(c,t,e,:)),squeeze(dec_values(c,t,e,:)),1); 
%         end
%     end
% end

%
% for z = 1:size(auc_pen,3)
%     auc_pen(:,:,z) = ndnanfilter(auc_pen(:,:,z),'rectwin',3);
% end

% return C and T with highest AUC
thresh = prctile(auc_pen(:),prctileThresh);
ind=find(auc_pen>=thresh,1,'first');
[cBest,tBest,encBest]=ind2sub(size(auc_pen),ind);
cBest = Cs(cBest);
encBest = types(encBest);



function [yProbs,yPreds,yTest,A,err,lambda,tBest,encBest,Cs,Ts,aucs] = doRegFun(X,Y,T,folds,iFold,lambda,normType,sessInds,trialType,prctileThresh,nestedCvGroup,percentCV,useKfold,Cs,tBest,encBest,doCalcPenalty)
% This does the classification.
% X = # trials x # features
% Y = # trials vector of responses
% folds = described above
% iFold = current fold number
% lambda = penalty parameter, if given
% normtype = L1 or L2
% sessInds = vector, same number of rows as X, identify the session for the
% row
%


% get train data for this fold
trainInds  = folds(iFold,:);
% if doPermute
%     randOrder = randperm(sum(trainInds));
%     yTmp = Y(trainInds);
%     yTmp = yTmp(randOrder);
%     Y(trainInds) = yTmp;
% end

sessions   = sessInds(trainInds);
% encType    = trialType(trainInds);
yTrainBool = Y(trainInds);
xTrain     = X(trainInds,:);

% if no lambda given, calculate lambda for this fold.
if doCalcPenalty
    if ~useKfold
        [~,subFolds] = createFolds(sessInds(trainInds)',nestedCvGroup(trainInds));
    else
        [~,subFolds] = createKfolds(sum(trainInds),percentCV);
    end
    [lambda,tBest,encBest,Cs,Ts,aucs] = calcPenalty(xTrain,yTrainBool,T,subFolds,sessions,trialType,normType,prctileThresh,Cs);
else    
    Ts = tBest;
    aucs = NaN;
end
    
if encBest == 3
    trialType(:) = 3;
end

% set parameters for either L1 or L2
if strcmpi(normType,'L1')
    liblin_param = ['-c ' sprintf('%f',lambda) ' -s 6 -q'];
elseif strcmpi(normType,'L2')
    liblin_param = ['-c ' sprintf('%f',lambda) ' -s 0 -q'];
end


% make sure the classes are balanced by doing a number of subsampling
% rounds. First see how many observations to remove to balance them
numToRemove = sum(yTrainBool) - sum(~yTrainBool)
maxSamp = 100;
% numToRemove = 0;
if numToRemove == 0
    maxSamp = 1;
end

% will store output from each round
preds = NaN(sum(~trainInds),maxSamp);
probs = NaN(sum(~trainInds),maxSamp);

for nSamp = 1:maxSamp
    
    % pick obs to remove
    toRemove = [];
    if numToRemove > 0
        toRemove = randsample(find(yTrainBool==1),abs(numToRemove));
    elseif numToRemove < 0
        toRemove = randsample(find(yTrainBool~=1),abs(numToRemove));
    end
    trainIndsSamp = setdiff(find(trainInds),toRemove);
    
    % remove from training set
    % training set x
    xTrain = X(trainIndsSamp,T==tBest & trialType==encBest);
    
    % training set y
    yTrain = double(Y(trainIndsSamp));
    yTrain(yTrain==0) = -1;
    
    % standardize training data
    [xTrain,m,s] = standardize(xTrain,sessInds(trainIndsSamp));
    
    % weight by class probabilities
    pos = mean(yTrain==1);
    neg = mean(yTrain==-1);
    mean_tmp = mean([pos neg]);
    pos = sprintf('%f',pos/mean_tmp);
    neg = sprintf('%f',neg/mean_tmp);
    param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
    
    % train model
    model = train(double(yTrain),sparse(xTrain),param);
    if nSamp == 1
        Ws = NaN(maxSamp,size(model.w,2));
    end
    Ws(nSamp,:) = model.w;
    
    % testing set
    xTest = X(~trainInds,T==tBest & trialType==encBest);
    yTest = double(Y(~trainInds));
    yTest(yTest==0) = -1;
    
    % standardize testing data
    xTest = standardize_test(xTest,sessInds(~trainInds),m,s);
    
    % predict
    [predSamp, acc, probSamp] = predict(double(yTest),sparse(xTest),model,'-b 1 -q');
    % same as p = glmval(model.w',xTest,'logit','constant','off')?
    
    preds(:,nSamp) = predSamp;
    probs(:,nSamp) = probSamp(:,1);
    
end

yProbs = mean(probs,2);
yPreds = mean(preds,2);
A      = mean(Ws,1);
err    = mean(preds == repmat(yTest,1,size(preds,2)),2);

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


function res = doRegFullData(X,Y,T,folds,normType,sessInds)
% train on all data. No test

% find lambda
res = [];
[lambda,bestT,Cs,Ts,aucs] = calcPenalty(X,Y,T,folds,sessInds,normType);
res.C = lambda;
res.T = bestT;
res.Cs = Cs;
res.Ts = Ts;
res.aucs = aucs;

% standardize
X = X(:,T==bestT);
[X,m,s] = standardize(X,sessInds);

% set parameters for either L1 or L2
if strcmpi(normType,'L1')
    liblin_param = ['-c ' sprintf('%f',lambda) ' -s 6'];
elseif strcmpi(normType,'L2')
    liblin_param = ['-c ' sprintf('%f',lambda) ' -s 0'];
end

% weight by class probabilities
Y = double(Y);
Y(Y==0) = -1;
pos = mean(Y==1);
neg = mean(Y==-1);
mean_tmp = mean([pos neg]);
pos = sprintf('%f',pos/mean_tmp);
neg = sprintf('%f',neg/mean_tmp);
param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];

% train model
model = train(double(Y),sparse(X),param);
res.A = model.w;


function [X,m,s] = standardize(X,sessInds)
% zscore the features

sessions = unique(sessInds);
m    = NaN(length(sessions),size(X,2));
s    = NaN(length(sessions),size(X,2));
for sess = 1:length(sessions)
    
    m(sess,:) = nanmean(X(sessInds==sessions(sess),:));
    s(sess,:) = nanstd(X(sessInds==sessions(sess),:));
    
    X(sessInds==sessions(sess),:) = bsxfun(@minus, X(sessInds==sessions(sess),:),...
        m(sess,:));
    
    X(sessInds==sessions(sess),:) = bsxfun(@rdivide, X(sessInds==sessions(sess),:),...
        s(sess,:));
end
% X(:,1) = 1;


function [X] = standardize_test(X,sessInds,m,s)
% zscore the features with given mean and standard deviation

sessions = unique(sessInds);
for sess = 1:length(sessions)
    
    X(sessInds==sessions(sess),:) = bsxfun(@minus, X(sessInds==sessions(sess),:),...
        m(sess,:));
    
    X(sessInds==sessions(sess),:) = bsxfun(@rdivide, X(sessInds==sessions(sess),:),...
        s(sess,:));
end
% X(:,1) = 1;

function [powerData,sessInds,timeLabel,tal] = loadWatrousFeatures(subj,params)

featureDir = fullfile('/home2/andrew.watrous/Results/YC1_Features',subj,'Features_4bands_all_elecs_11_16_2015');
sessions   = dir(fullfile(featureDir,['Session_*']));


subjDataPow   = cell(1,4);
subjDataPhase = cell(1,4);
sessInds = [];
for s = 1:length(sessions)
    fname = fullfile(featureDir,sessions(s).name,'Subject_features.mat');
    sessData = load(fname);
    [tal,elecs] = filterTalByRegion(sessData.features.elecs,params.region);
    sessInds = cat(1,sessInds,ones(size(sessData.features.pow.band{1},3),1)*s);
    for band = 1:4
        subjDataPow{band} = cat(3,subjDataPow{band},sessData.features.pow.band{band}(elecs,:,:));
        subjDataPhase{band} = cat(3,subjDataPhase{band},sessData.features.phase.band{band}(elecs,:,:));
    end
end

nTimes = size(params.timeBins,1);

dataPow = NaN(length(sessInds),4,nTimes,size(subjDataPow{1},1));
dataPhase = NaN(length(sessInds),4*2,nTimes,size(subjDataPhase{1},1));
timeLabelPow = NaN(1,nTimes,size(dataPow,3),size(dataPow,4));
timeLabelPhase = NaN(1,nTimes,size(dataPhase,3),size(dataPhase,4));


for f = 1:4
    bandPow = subjDataPow{f};
    bandPow = permute(bandPow,[3 2 1]);
    
    bandPhase = subjDataPhase{f};
    bandPhase = permute(bandPhase,[3 2 1]);
    
    for t = 1:nTimes
        tStart = params.timeBins(t,1);
        tEnd = params.timeBins(t,2);
        timeBandPhase = bandPhase(:,tStart:tEnd,:);
        
        % can't specific a dimenstion in circmean, so have to loop?
        for ev = 1:size(timeBandPhase,1);
            for elec = 1:size(timeBandPhase,3)
                theta = circmean(timeBandPhase(ev,:,elec));
                dataPhase(ev,f*2-1,t,elec) = sin(theta);
                dataPhase(ev,f*2,t,elec) = cos(theta);
            end
        end
        
        timeBandPow = nanmean(bandPow(:,tStart:tEnd,:),2);
        dataPow(:,f,t,:) = timeBandPow;
        timeLabelPow(:,:,t,:) = t;
        timeLabelPhase(:,:,t,:) = t;
    end
end

if params.usePhase==1
    powerData = dataPhase;
    timeLabel = timeLabelPhase;
elseif params.usePhase==0
    powerData = dataPow;
    timeLabel = timeLabelPow;
elseif params.usePhase == 2
    powerData = cat(2,dataPow,dataPhase);
    timeLabel = cat(2,timeLabelPow,timeLabelPhase);
end


% tal = sessData.features.elecs(elecs);
function [powerData,sessInds] = loadAllPower(tal,subj,events,timeBins,eventsToUse,params)


nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
nElecs = length(tal);


% when loading power, use either original power or power with effect of
% trial number removed.
powField  = 'pow';
stdField  = 'stdBasePow';
meanField = 'meanBasePow';
if params.useCorrectedPower
    powField  = 'powCorr';
    stdField  = 'stdBasePowCorr';
    meanField = 'meanBasePowCorr';
end
for e = 1:nElecs
    elecNum = tal(e).channel;
    
    basePath  = params.powerPath;
    subjPath  = fullfile(basePath,subj);
    sessions  = unique([events.session]);
    subjPow   = [];
    subjPhase = [];
    sessInds  = [];
    
    for s = 1:length(sessions)
        fname = fullfile(subjPath,'RAM_YC1_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
        sessPow = load(fname);
        
        if e == 1;
            nFreqs = size(sessPow.sessOutput.pow,1);
            powerData = NaN(nFreqs,nTimes,nEvents,nElecs);
            phaseData = NaN(nFreqs*2,nTimes,nEvents,nElecs);
        end
        
        % under zscore
        zPow = sessPow.sessOutput.(powField);
        nT = size(sessPow.sessOutput.(powField),2);
        nE = size(sessPow.sessOutput.(powField),3);
        pow = zPow .* repmat(sessPow.sessOutput.(stdField),[1 nT nE]) + repmat(sessPow.sessOutput.(meanField),[1 nT nE]);
        %         pow = zPow;
        subjPow = cat(3,subjPow,pow);
        subjPhase = cat(3,subjPhase,sessPow.sessOutput.phase);
        sessInds = cat(1,sessInds,ones(nE,1)*s);
    end
    
    if length(eventsToUse) ~= size(subjPow,3)
        fprintf('Number of events does not match size of power matrix for %s!.\n',subj)
        return
    end
    subjPow    = subjPow(:,:,eventsToUse);
    subjPhase  = subjPhase(:,:,eventsToUse);
    sessInds   = sessInds(eventsToUse);
    
    % average times
    for t = 1:nTimes
        tStart = params.timeBins(t,1);
        tEnd = params.timeBins(t,2);
        tInds = tStart:tEnd;
        powerData(:,t,:,e) = nanmean(subjPow(:,tInds,:),2); 
                
        % can't specific a dimenstion in circmean, so have to loop?
        for ev = 1:size(phaseData,3);
            for f = 1:nFreqs
                theta = circmean(subjPhase(f,tInds,ev));
                phaseData(f*2-1,t,ev,e) = sin(theta);
                phaseData(f*2,t,ev,e) = cos(theta);
            end
        end 
    end    
end

if params.usePhase==1
    powerData = phaseData;    
elseif params.usePhase == 2
    powerData = cat(1,powerData,phaseData);    
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














