function [perf,AUC,subject,params] = YC1_runMulti_subj_combinedFeatures(subj,params,saveDir)
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
%      AUC: area under the ROC curve
%     perf: percent classifier accuracy
%  subject: current subject
%

perf    = [];
subject = [];
AUC     = [];


% do we overwrite?
fname = fullfile(saveDir,[subj '_lasso.mat']);
if exist(fname,'file') && params.saveOutput && ~params.overwrite
    fprintf('Lasso file already exists for %s.\n',subj)
    return
end

try
    
    % load subject electrode locations and filter to specific regions if
    % desired.
    tal = getBipolarSubjElecs(subj,1,1,params.excludeEpiElecs);
    tal = filterTalByRegion(tal,params.region);
    if isempty(tal)
        fprintf('No %s electrode for %s.\n',params.region,subj)
        return
    end
    
    % load power parameters
    powParams = load(fullfile(params.powerPath,'params_RAM_YC1.mat'));
    
    % Set up time bins
    tEnds     = (powParams.params.pow.timeWin:powParams.params.pow.timeStep:powParams.params.eeg.durationMS)+powParams.params.eeg.offsetMS;
    tStarts   = tEnds - powParams.params.pow.timeWin+1;
    powParams.timeBins = [tStarts' tEnds'];
    
    % get parameters
    freqBins      = params.freqBins;
    timeBins      = params.timeBins;
    modelEachTime = params.modelEachTime;
    doBinary      = params.doBinary;
    saveOutput    = params.saveOutput;
    doPermute     = params.doPermute;
    normType      = params.normType;
    useKfold      = params.useKfold;
    cvField       = params.cvField;   
    nestedCvField = params.nestedCvField; 
    prctileThresh = params.auc_prctileThresh;
    
    % load events
    events = get_sub_events('RAM_YC1',subj);
    
    % add the test error to the learning trials
    events    = addExtraFields(events);
    session   = [events.session];    
    
    % filter to events of interest
    eventsToUse = params.eventFilter(events);
    if strcmpi(cvField,'session') && length(unique([events.session])) < 2
        fprintf('Not enough sessions for %s using session cvField.\n',subj)
        return
    end
    if sum(eventsToUse) < 10
        fprintf('Not enough events for %s.\n',subj)
        return
    end        
    encPeriod = [events(eventsToUse).withinItemCount];
    
    % load power for all electrodes
    powDir = fullfile(saveDir,'power');
    powFile = fullfile(powDir,[subj '_binnedPower.mat']);
    if ~exist(powDir,'dir')
        mkdir(powDir)
    end
    if params.loadPower
        powerData = load(powFile);        
        timeLabel = powerData.timeLabel;
        encLabel = powerData.encLabel;
        sessInds = powerData.sessInds;
        powerData = powerData.powerData;
    else
        
        % load the power
        [powerData,sessInds] = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params);
        powerData = permute(powerData,[3 1 2 4]);
        
        % matrix to keep track of encoding trial 1 or 2 feature        
        firstEnc  = powerData(encPeriod==1,:,:,:);
        secondEnc = powerData(encPeriod==2,:,:,:);
        powerData = cat(5,firstEnc,secondEnc);
        encLabel = NaN(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5));
        for t = 1:2
            encLabel(:,:,:,:,t) = t;
        end        
        
        % create a matrix to keep track of the time bin features
        timeLabel = NaN(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5));        
        for t = 1:size(timeLabel,3)
            timeLabel(:,:,t,:) = t;
        end
                
        if params.savePower
            powFile = fullfile(powDir,[subj '_binnedPower.mat']);
            save(powFile,'powerData','params','sessInds','timeLabel','encLabel')
        end
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
%     if doPermute
%         randOrder = randperm(size(trials,1));        
%         if ~isfield(params,'encPeriod') || strcmpi(params.encPeriod,'both')
%             randOrder = [randOrder;randOrder];
%             randOrder = randOrder(:);
%             Y = Y(randOrder);
%         else
%             Y = Y(randOrder);
%         end        
%     end
          
    % Unlike previous versions, if modelEachTime is true, this will not
    % create nTimes number of models. Instead, the model for fold each will
    % be chosen based on inner fold cross validation of both penalty and
    % time point.
    res = [];
    X   = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs*nItemObs);    
%     m = repmat(nanmean(X,1),size(X,1),1); 
%     X = X - m;
    
    T   = reshape(squeeze(timeLabel),size(timeLabel,1),nFreqs*nTimes*nElecs*nItemObs);    
%     X   = standardize(X,sessInds);
    
    % will also do CV for encoding period 1 or 2
    trialType = reshape(squeeze(encLabel),size(encLabel,1),nFreqs*nTimes*nElecs*nItemObs);
    
    % see if we are precomputing lambda based on all the data. Again,
    % we shouldn't really do this
    lambda = [];
    if params.crossValStrictness == 0
        if ~isempty(params.lambda)
            lambda = params.lambda;
        else
            fprintf('Subject %s: Computing optimal c using %s regularization.\n',subj,normType)
            lambda = calcPenalty(X,Y,T,folds,sessInds,normType);
        end
    end
    
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
            res.aucs{iFold}] = doRegFun(X,Y,T,folds,iFold,lambda,normType,sessInds,trialType,prctileThresh,nestedCvGroup,params.percentCV,useKfold,doPermute);
        fprintf('Subject %s: Fold %d of %d, %.3f running percent correct.\n',subj,iFold,nFolds,mean(vertcat(res.err{:})))
        
    end
    
    yTest = vertcat(res.yTest{:});
    yProb = vertcat(res.yProb{:});
    [~,~,~,res.AUC] = perfcurve(yTest,yProb,1);
    AUC = res.AUC;
    perf = mean(vertcat(res.err{:}));
    
    % Finally, create a model based on the full dataset for the patient.
%     resFull = doRegFullData(X,Y,T,folds,normType,sessInds,trialType,prctileThresh);
        

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


function [cBest,tBest,encBest,Cs,Ts,auc_pen] = calcPenalty(X,Y,T,folds,sessInds,trialType,normType,prctileThresh)
% Calculate the optimal penalty parameter.
%
% Returns: penalty (the optimal C)
%          C (vector C values tried)
%          auc_pen (the auc for each penalty paramter)

Y = double(Y);
Y(Y==0) = -1;
Cs      = logspace(log10(1e-1),log10(1e4),22);
Ts      = 1:length(unique(T));
nCV     = size(folds,1);
types   = unique(trialType);
auc_pen = NaN(length(Cs),length(Ts),length(types)+1);
% auc_pen = NaN(1,length(Cs)*length(Ts)*(length(types)+1));

maxInd = length(Cs)*length(Ts)*(length(types)+1);

parfor ind = 1:maxInd
    
    [thisC,thisT,thisEnc] = ind2sub([length(Cs),length(Ts),length(types)+1],ind);
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
    dec_values = NaN(nCV,1);
    labels = NaN(nCV,1);
    
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
        
        % predict with model
        [pred, acc, dec] = predict(double(yTest),sparse(xTest),model,'-b 1 -q');
        
        % keep of track of output over all folds
%         if model.Label(1) < 0;
%             dec = dec * -1;
%         end
        dec_values(cv) = dec(:,1);
        labels(cv) = yTest;               
    end
    [~,~,~,auc_pen(ind)] = perfcurve(labels,dec_values,1);
end

for z = 1:size(auc_pen,3)
    auc_pen(:,:,z) = ndnanfilter(auc_pen(:,:,z),'rectwin',3);
end

% return C and T with highest AUC
thresh = prctile(auc_pen(:),prctileThresh);
ind=find(auc_pen>=thresh,1,'first');
[cBest,tBest,encBest]=ind2sub(size(auc_pen),ind);
% tBest = 11
 
% tBest = 11; 
% encBest = 3; 
% cBest = 6;
% r = randperm(length(Cs));
% cBest = r(1);
% cBest = 10
% r = randperm(length(Ts));
% tBest = r(1);
% tBest = 11
% r = randperm(3);
% encBest = r(1);
% enBest = 2
cBest = Cs(cBest);


function [yProbs,yPreds,yTest,A,err,lambda,tBest,encBest,Cs,Ts,aucs] = doRegFun(X,Y,T,folds,iFold,lambda,normType,sessInds,trialType,prctileThresh,nestedCvGroup,percentCV,useKfold,doPermute)
% This does the classification.
%
% Inputs:
%
%             X: # trials x # features
%             Y: # trials vector of responses
%             T: # of features long vector, indicating the timepoint the
%                feature is from
%         folds: described above
%         iFold: current fold number
%        lambda: penalty parameter, if given
%      normtype: L1 or L2
%      sessInds: vector, same number of rows as X, identify the session
%                for the row
%        events: events struture, one entry for each obcervation
%  nestedCvFied: grouping field of the eventus structure from which to
%                create the nested cross val folds
% modelEachTime:
%
%
% Outputs:
%
% yProbs
% yTest
% A
% err
% lambda
% tBest
% Cs
% Ts
% aucs



% get train data for this fold
trainInds  = folds(iFold,:);
if doPermute               
    randOrder = randperm(sum(trainInds));
    yTmp = Y(trainInds);
    yTmp = yTmp(randOrder);
    Y(trainInds) = yTmp;
end

sessions   = sessInds(trainInds);
% encType    = trialType(trainInds);
yTrainBool = Y(trainInds);
xTrain     = X(trainInds,:);

% if no lambda given, calculate lambda for this fold.
if isempty(lambda)    
    % trainInds  = folds(iFold,setdiff(1:size(folds,2),toRemove));
    % subFolds   = folds(setdiff(1:size(folds,1),iFold),trainInds);    
    
    if ~useKfold
        [~,subFolds] = createFolds(sessInds(trainInds)',nestedCvGroup(trainInds));
    else
        [~,subFolds] = createKfolds(sum(trainInds),percentCV);
    end
    [lambda,tBest,encBest,Cs,Ts,aucs] = calcPenalty(xTrain,yTrainBool,T,subFolds,sessions,trialType,normType,prctileThresh);         
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

% use only the enc period chosen by the cross validation
% X          = X(trialType==encBest,:);
% Y          = Y(trialType==encBest);
% trainInds  = trainInds(trialType==encBest);
% yTrainBool = Y(trainInds);


% make sure the classes are balanced by doing a number of subsampling
% rounds. First see how many observations to remove to balance them
numToRemove = sum(yTrainBool) - sum(~yTrainBool);
maxSamp = 100;
numToRemove = 0;
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
    %     mTrain = nanmean(xTrain,1);
    %     m = repmat(mTrain,size(xTrain,1),1);
    %     xTrain = xTrain - m;
    
    % training set y
    yTrain = double(Y(trainIndsSamp));
    yTrain(yTrain==0) = -1;
    
    % standardize training data
    % xTrain = [ones(size(xTrain,1),1) xTrain];
        [xTrain,m,s] = standardize(xTrain,sessInds(trainIndsSamp));
%     [xTrain,xMin,xMax] = normalize(xTrain,sessInds(trainIndsSamp));
    
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
    % xTest = [ones(size(xTest,1),1) xTest];
%     m = repmat(mTrain,size(xTest,1),1);
%     xTest = xTest - m;
    xTest = standardize_test(xTest,sessInds(~trainInds),m,s);
%      xTest = normalize_test(xTest,sessInds(~trainInds),xMin,xMax);
    
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

function res = doRegFullData(X,Y,T,folds,normType,sessInds,trialType,prctileThresh)
% train on all data. No test

% find lambda
res = [];
[lambda,bestT,encBest,Cs,Ts,aucs] = calcPenalty(X,Y,T,folds,sessInds,trialType,normType,prctileThresh); 
res.C = lambda;
res.T = bestT;
res.Cs = Cs;
res.Ts = Ts;
res.encBest = encBest;
res.aucs = aucs;

% standardize 
X = X(:,T==bestT&trialType==encBest);
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

function [X,xMin,xMax] = normalize(X,sessInds)
% normalize the features

sessions = unique(sessInds);
xMin    = NaN(length(sessions),size(X,2));
xMax    = NaN(length(sessions),size(X,2));

for sess = 1:length(sessions)
   
    xMin(sess,:) = min(X(sessInds==sessions(sess),:));
    xMax(sess,:) = max(X(sessInds==sessions(sess),:));    
    
    
    
    sessData = X(sessInds==sessions(sess),:);
%     sessData = bsxfun(@minus,sessData,xMin(sess));
%     sessData = bsxfun(@times,sessData,2./xMax(sess));    
    sessDataNorm = NaN(size(sessData));
    for col = 1:size(sessData,2);    
        sessDataNorm(:,col) = ((sessData(:,col) - xMin(sess,col))./(xMax(sess,col)-xMin(sess,col)) - .5)*2;
%         sessDataNorm(:,col) = 2*mat2gray(sessData(:,col))-1;        
    end
    X(sessInds==sessions(sess),:) = sessDataNorm;    
end
% X(:,1) = 1;

function [X] = normalize_test(X,sessInds,xMin,xMax)

sessions = unique(sessInds);
for sess = 1:length(sessions)
 
    sessData = X(sessInds==sessions(sess),:);
    sessDataNorm = NaN(size(sessData));
    for col = 1:size(sessData,2);    
        sessDataNorm(:,col) = ((sessData(:,col) - xMin(sess,col))./(xMax(sess,col)-xMin(sess,col)) - .5)*2;
%         sessDataNorm(:,col) = 2*mat2gray(sessData(:,col),[xMin(col) xMax(col)])-1;        
    end
    X(sessInds==sessions(sess),:) = sessDataNorm;    
end


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



function [powerData,sessInds] = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params)

nFreqs = size(freqBins,1);
nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
nElecs = length(tal);
powerData = NaN(nFreqs,nTimes,nEvents,nElecs);

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
    sessions = unique([events.session]);
    subjPow  = [];
    sessInds = [];
    
    for s = 1:length(sessions)
        fname = fullfile(subjPath,'RAM_YC1_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
        sessPow = load(fname);
        
        % under zscore
        zPow = sessPow.sessOutput.(powField);
        nT = size(sessPow.sessOutput.(powField),2);
        nE = size(sessPow.sessOutput.(powField),3);
        pow = zPow .* repmat(sessPow.sessOutput.(stdField),[1 nT nE]) + repmat(sessPow.sessOutput.(meanField),[1 nT nE]);    
%         pow = zPow;
        subjPow = cat(3,subjPow,pow);
        sessInds = cat(1,sessInds,ones(nE,1)*s);
    end
    
    if length(eventsToUse) ~= size(subjPow,3)
        fprintf('Number of events does not match size of power matrix for %s!.\n',subj)
        return
    end
    subjPow  = subjPow(:,:,eventsToUse);
    sessInds = sessInds(eventsToUse);    
        
    
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
        tInds = powParams.timeBins(:,1) >= timeBins(t,1) & powParams.timeBins(:,2) <= timeBins(t,2);
        tmpPower(:,t,:) = nanmean(subjPow(:,tInds,:),2);
    end
    powerData(:,:,:,e) = tmpPower;
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











