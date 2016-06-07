function [perf,AUC,subject,params,res] = YC1_watrous_phase(subj,params,saveDir)
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
res     = [];

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
    [powMat,~,tal,powParams] = loadAllPower(tal,subj,events,params.timeBins,eventsToUse,params);
    nElecs = size(powMat,2);
    nEncs  = size(powMat,3);
    nTimes = size(powMat,4);
    nFreqs = size(powMat,5);
    

%     % create matrix to keep track of time bins
%     timeMat = NaN(1,nElecs,nEncs,nTimes,nFreqs);
%     for t = 1:size(timeMat,4)
%         timeMat(:,:,:,t,:) = t;
%     end
%     
%     % create matrix to keep track of enc periods
%     encMat = NaN(1,nElecs,nEncs,nTimes,nFreqs);
%     for e = 1:size(encMat,3)
%         timeMat(:,:,:,3,:) = e;
%     end
    
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
        [yProbFold,yTestFold] = run_fold_with_nested_folds(powMat,Y,sessions,folds,thisFold,sel,params,powParams);        
        yProbs = [yProbs,yProbFold];
        yTests = [yTests,yTestFold];
        if ismember([-1,1],yTests)
            [~,~,~,aucTmp] = perfcurve(yTests,yProbs,1);
            fprintf('%s: running AUC = %.3f\n',subj,aucTmp)
        end
       
    end
    
    auc = aucTmp;
    save(fname,'yProbs','yTests','auc')
    
%     
%     % if we are given one penality parameter and one time bin, and we are
%     % not cross validating the encoding period, then we don't need an inner
%     % CV loop
%     tBest = [];
%     encBest = [];
%     C = [];
%     G = [];
%     doCalcPenalty = 1;
%     if strcmp(normType,'svm')
%         if length(Cs) == 1 && length(Gs) == 1 && size(params.timeBins,1) == 1 && ~strcmp(params.encPeriod,'both')
%             doCalcPenalty = 0;
%             C = Cs;
%             G = Gs;
%             tBest = 1;
%             if strcmp(params.encPeriod,'first')
%                 encBest = 1;
%             elseif strcmp(params.encPeriod,'second')
%                 encBest = 2;
%             elseif strcmp(params.encPeriod,'combined')
%                 encBest = 3;
%             end
%         end
%     else
%         if length(Cs) == 1 && size(params.timeBins,1) == 1 && ~strcmp(params.encPeriod,'both')
%             doCalcPenalty = 0;
%             C = Cs;
%             tBest = 1;
%             if strcmp(params.encPeriod,'first')
%                 encBest = 1;
%             elseif strcmp(params.encPeriod,'second')
%                 encBest = 2;
%             elseif strcmp(params.encPeriod,'combined')
%                 encBest = 3;
%             end
%         end
%     end
%     
%     
% 
%     % SOME CHECKS TO ENSURE ENOUGH EVENTS AND SESSIONS
%     %---------------------------------------------------------------------%
%     eventsToUse = params.eventFilter(events);
%     if sum(eventsToUse) < 10
%         fprintf('Not enough events for %s.\n',subj)
%         return
%     end
%     encPeriod = [events(eventsToUse).withinItemCount];
%     
%     if strcmpi(nestedCvField,'session')
%         if (length(unique(session)) <= 2 && doCalcPenalty) || (~doCalcPenalty && length(unique(session)) == 1)
%             fprintf('Not enough sessions for session nested crossval for %s.\n',subj)
%             return
%         end
%     end
%     if strcmpi(cvField,'session')
%         if length(unique(session)) < 2
%             fprintf('Not enough sessions for session outer crossval for %s.\n',subj)
%             return
%         end
%     end
%     %---------------------------------------------------------------------%
%     
%     % load power or phase
%     if params.useWatrous == 1
%         [powerData,sessInds,~,tal] = loadWatrousFeatures(subj,params);
%     else
%         tal = getBipolarSubjElecs(subj,1,1,params.excludeEpiElecs);
%         [powerData,sessInds,tal] = loadAllPower(tal,subj,events,params.timeBins,eventsToUse,params);
%         powerData = permute(powerData,[3 1 2 4]);
%     end
%     
%     if isempty(tal)
%         fprintf('no mtl elecs for %s.\n',subj)
%         return
%     end
%     
%     % make sure the power values correspond to the events I expect
%     if sum(eventsToUse) ~= size(powerData,1)
%         fprintf('EVENTS MISMATCH FOR %s.\n',subj)
%         return
%     end
%     
%     % matrix to keep track of encoding trial 1 or 2 feature
%     firstEnc  = powerData(encPeriod==1,:,:,:);
%     secondEnc = powerData(encPeriod==2,:,:,:);
%     
%     if strcmp(params.encPeriod,'first')
%         powerData = firstEnc;
%         encLabel  = ones(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5));
%     elseif strcmp(params.encPeriod,'second')
%         powerData = secondEnc;
%         encLabel  = ones(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5))+1;
%     elseif strcmp(params.encPeriod,'both')
%         powerData = cat(5,firstEnc,secondEnc);
%         encLabel = NaN(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5));
%         for t = 1:2
%             encLabel(:,:,:,:,t) = t;
%         end
%     elseif strcmp(params.encPeriod,'combined')
%         powerData = cat(5,firstEnc,secondEnc);
%         encLabel = ones(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5))+2;
%     end
%     
%     % create a matrix to keep track of the time bin features
%     timeLabel = NaN(1,size(powerData,2),size(powerData,3),size(powerData,4),size(powerData,5));
%     for t = 1:size(timeLabel,3)
%         timeLabel(:,:,t,:,:) = t;
%     end
%     
%     % size of feature matrix
%     nElecs = size(powerData,4);
%     nTimes = size(powerData,3);
%     nFreqs = size(powerData,2);
%     nItemObs = size(powerData,5);
%     
%     % response data
%     sessInds = sessInds(encPeriod==1);
%     if isfield(params,'diffRegress') && params.diffRegress == 1
%         Y = [events(eventsToUse).residual]';
%         Y = Y(encPeriod==1);
%         Y = Y > 0;
%     else
%         Y = [events(eventsToUse).testError]';
%         Y = Y(encPeriod==1);
%         Y = Y < median(Y);
% %         Y = Y < .3;
%     end
%     
%     
%     % DON'T FORGET TO TAEK THIS OUT!!!!!!
% %     Y = Y < .3;
% %     if doBinary
% %         Y  = Y < median(Y);
% %     end
%     
%     % determine the cross validation folds.
%     session = session(eventsToUse);
%     cvGroup = [events(eventsToUse).(cvField)];
%     nestedCvGroup = [events(eventsToUse).(nestedCvField)];
%     nestedCvGroup = nestedCvGroup(encPeriod==1);
%     
%     if ~useKfold
%         [nFolds,folds,trials] = createFolds(session(encPeriod==1),cvGroup(encPeriod==1));
%     else
%         [nFolds,folds,trials] = createKfolds(sum(encPeriod==1),params.percentCV);
%     end
%     
%     % permute the responses if desired
%     
% %     if doPermute
% %         randOrder = randperm(size(trials,1));
% %         if ~isfield(params,'encPeriod') || strcmpi(params.encPeriod,'both')
% %             randOrder = [randOrder;randOrder];
% %             randOrder = randOrder(:);
% %             Y = Y(randOrder);
% %         else
% %             Y = Y(randOrder);
% %         end
% %     end
%     
%     
%     % reshape into obs x features
%     res = [];
%     X         = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs*nItemObs);
%     X         = standardize(X,sessInds);
%     T         = reshape(squeeze(timeLabel),size(timeLabel,1),nFreqs*nTimes*nElecs*nItemObs);
%     trialType = reshape(squeeze(encLabel),size(encLabel,1),nFreqs*nTimes*nElecs*nItemObs);
% %     if params.randX
% %         X = randn(size(X));
% %     end
%     
% %     fakeX = NaN(size(X));
% %     for row = 1:size(X,1)
% %         for c = 1:size(X,2)
% %             
% %             r = randi(size(X,1));
% %             fakeX(row,c) = X(r,c);
% %         end
% %     end
% %     X = fakeX;
%  
%     % run for each fold
%     [res.yProb,res.yPred,res.yTest,res.A,res.intercept,res.err,...
%         res.C,res.G,res.tBest,res.encBest,res.Cs,res.Gs,res.Ts,res.aucs,res.fold_auc] = deal(cell(nFolds,1));
%     for iFold = 1:nFolds
%         fprintf('Subject %s: Fold %d of %d.\n',subj,iFold,nFolds)
%         [res.yProb{iFold},...
%             res.yPred{iFold},...
%             res.yTest{iFold},...
%             res.A{iFold},...
%             res.err{iFold},...
%             res.C{iFold},...
%             res.G{iFold},...
%             res.tBest{iFold},...
%             res.encBest{iFold},...
%             res.Cs{iFold},...
%             res.Gs{iFold},...
%             res.Ts{iFold},...
%             res.aucs{iFold}] = doRegFun(X,Y,T,folds,iFold,C,G,normType,sessInds,trialType,prctileThresh,nestedCvGroup,params.percentCV,useKfold,Cs,Gs,tBest,encBest,doCalcPenalty,doPermute);
%         if strcmpi(cvField,'session')                        
%             [~,~,~,res.fold_auc{iFold}] = perfcurve(res.yTest{iFold},res.yProb{iFold},1);
%         end
%         fprintf('Subject %s: Fold %d of %d, %.3f running percent correct.\n',subj,iFold,nFolds,mean(vertcat(res.err{:})))        
%     end
%     
%     yTest = vertcat(res.yTest{:});
%     yProb = vertcat(res.yProb{:});
%     [~,~,~,res.AUC] = perfcurve(yTest,yProb,1);
%     AUC = res.AUC;
%     perf = mean(vertcat(res.err{:}));
%     
%     % Finally, create a model based on the full dataset for the patient. If
%     % modelEachTime, do it for the optimal time
%     %     bestTime = find(AUC == max(AUC),1,'first');
%     %     X = reshape(squeeze(powerData(:,:,bestTime,:,:)),size(powerData,1),nFreqs*nElecs*nItemObs);
%     %     resFull = doRegFullData(X,Y,T,folds,normType,sessInds);
%     
%     subject       = subj;
%     params.C = C;
%     if saveOutput
%         objLocs = vertcat(events(eventsToUse).objLocs);
%         save(fname,'res','Y','objLocs','params','perf','tal','AUC');
%     end
catch e
    fname = fullfile(saveDir,[subj '_class_error.mat']);
    save(fname,'e')
end

function [yProb,yTest] = run_fold_with_nested_folds(X,Y,sessions,folds,thisFold,sel,params,powParams)

% convert Y into double for liblinear
Y = double(Y);
Y(Y==0) = -1;

% unique times and enc periods
Ts = 1:size(X,4);
Es = params.encPeriods;
Cs = params.Cs;

% get training data
trainMask = folds(thisFold,:);
[~,subFolds] = createFolds(sessions(trainMask),sel(trainMask));
xTrain    = X(trainMask,:,:,:,:);
yTrain    = Y(trainMask);
sessTrain = sessions(trainMask);

% run nested fold to find optimal enc, time bin, and penalty
aucs = run_nested_fold(xTrain,yTrain,sessTrain,subFolds,params,powParams);
ind=find(aucs==max(aucs(:)),1,'first');
[eBestInd,tBestInd,cBestInd]=ind2sub(size(aucs),ind);
thisE = Es(eBestInd);
thisT = Ts(tBestInd);
thisC = Cs(cBestInd);
fprintf('Max AUC = %.3f, best enc = %d, best time = %d, best C = %.4f\n',max(aucs(:)),thisE,thisT,thisC)

% normalize full training set        
[xTrain,T,E,m,s,sessMap] = normPow_andMeanFreqs(xTrain,sessTrain,params,powParams.params.pow.freqs,[],[],[]);

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
param = liblin_param;
model = trainFun(double(yTrain),sparse(xTrain),param);

% test data
xTest    = X(~trainMask,:,:,:,:);
yTest    = Y(~trainMask);
sessTest = sessions(~trainMask);

% normalize test data and mean into freqs of interest
xTest = normPow_andMeanFreqs(xTest,sessTest,params,powParams.params.pow.freqs,m,s,sessMap);

% test on held out
xTest = xTest(:,featureMask);
[pred, acc, dec] = testFun(double(yTest),sparse(xTest),model,'-b 1 -q');

% probability of positive class
yProb = dec(:,model.Label==1);

function aucs = run_nested_fold(X,Y,sessions,subFolds,params,powParams)

% convert Y into double for liblinear
Y = double(Y);
Y(Y==0) = -1;

% unique times and enc periods
Ts = 1:size(X,4);
Es = params.encPeriods;
Cs = params.Cs;

% will hold auc for each combination of parametres
aucs = NaN(length(Es),length(Ts),length(Cs));

% max loop number for parfor
maxInd = length(Cs)*length(Ts)*length(Es);

parfor ind = 1:maxInd
    
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
    
    dec_values = [];
    labels = [];
    for thisFold = 1:size(subFolds,1)
        
        % masks for train observations
        obsMask  = subFolds(thisFold,:);        
        
        % training data for fold
        xTrain    = X(obsMask,:,:,:,:);
        yTrain    = Y(obsMask);
        sessTrain = sessions(obsMask);
        
        % normalize training data and mean into frequencies of interest         
        [xTrain,T,E,m,s,sessMap] = normPow_andMeanFreqs(xTrain,sessTrain,params,powParams.params.pow.freqs,[],[],[]);
        
        % mask encoding period and time bin
        featureMask = E == thisE & T == thisT;
        xTrain = xTrain(:,featureMask);
        
        % weight by percentage of positive and negative classes
        % pos = mean(yTrain==1);
        % neg = mean(yTrain==-1);
        % mean_tmp = mean([pos neg]);
        % pos = sprintf('%f',pos/mean_tmp);
        % neg = sprintf('%f',neg/mean_tmp);
        % param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
                
        % create model
        param = liblin_param;
        model = trainFun(double(yTrain),sparse(xTrain),param);
        
        % test data
        xTest    = X(~obsMask,:,:,:,:);
        yTest    = Y(~obsMask);
        sessTest = sessions(~obsMask);
        
        % normalize test data and mean into freqs of interest
        xTest = normPow_andMeanFreqs(xTest,sessTest,params,powParams.params.pow.freqs,m,s,sessMap);        
        
        % test on held out
        xTest = xTest(:,featureMask);               
        [pred, acc, dec] = testFun(double(yTest),sparse(xTest),model,'-b 1 -q');
        
        dec_values = [dec_values; dec(:,model.Label==1)];
        labels = [labels;yTest];        
    end  
    [~,~,~,aucs(ind)] = perfcurve(labels,dec_values,1);
end


function [cBest,gBest,tBest,encBest,Cs,Gs,Ts,auc_pen] = calcPenalty(X,Y,T,folds,sessInds,trialType,normType,prctileThresh,Cs,Gs)
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
    elseif strcmpi(normType,'svm')
        Cs = 10.^(-2:5);
        Cs = 2.^(-5:15);
    end
end
if isempty(Gs) && strcmpi(normType,'svm')
    Gs = 10.^(-5:1);
    Gs = 2.^(-15:3);
end
Ts      = 1:length(unique(T));
nCV     = size(folds,1);
types   = unique(trialType);

if strcmpi(normType,'svm')
    auc_pen = NaN(length(Cs),length(Ts),length(types),length(Gs));
    maxInd = length(Cs)*length(Ts)*length(types)*length(Gs);
else
    auc_pen = NaN(length(Cs),length(Ts),length(types));
    maxInd = length(Cs)*length(Ts)*length(types);
end
 
parfor ind = 1:maxInd
   
    if strcmpi(normType,'svm')
        [thisC,thisT,thisEnc,thisG] = ind2sub([length(Cs),length(Ts),length(types),length(Gs)],ind);
        thisKern = Gs(thisG);
    else
        [thisC,thisT,thisEnc] = ind2sub([length(Cs),length(Ts),length(types)],ind);
    end
    thisPen = Cs(thisC);
    
    % set parameters for either L1 or L2 with the current c
    if strcmpi(normType,'L1')
        liblin_param = ['-c ' sprintf('%f',thisPen) ' -s 6 -q'];
        trainFun = @train;
        testFun  = @predict;
    elseif strcmpi(normType,'L2')
        liblin_param = ['-c ' sprintf('%f',thisPen) ' -s 0 -q'];
        trainFun = @train;
        testFun  = @predict;
    elseif strcmpi(normType,'svm')
        liblin_param = ['-t 2 -g ' sprintf('%f',thisKern) ' -c ' sprintf('%f',thisPen) ' -b 0 -q'];
        trainFun = @svmtrain;
        testFun  = @svmpredict;
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
        
        % make sure the classes are balanced by doing a number of subsampling
        % rounds. First see how many observations to remove to balance them
%         keyboard
        yTrainBool = Y(inds);
        numToRemove = sum(yTrainBool==1) - sum(yTrainBool==-1);
        numToRemove = 0;
        toRemove = [];
%         if numToRemove~=0;keyboard
%         end
        if numToRemove > 0
            toRemove = randsample(unique([find(~inds) find(Y==1)']),abs(numToRemove));
        elseif numToRemove < 0
            toRemove = randsample(unique([find(~inds) find(Y~=1)']),abs(numToRemove));
        end
        trainIndsSamp = setdiff(find(inds),toRemove);
        
        % train data for this cv
        xTrain = X(trainIndsSamp,tInds&encInds);
        [xTrain,m,s] = standardize(xTrain,ones(size(xTrain,1),1));
        yTrain = Y(trainIndsSamp);
%         [xTrain,m,s] = standardize(xTrain,sessInds(trainIndsSamp));        
        
        % test data for this cv
        xTest =  X(~inds,tInds&encInds);
        %         xTest = [ones(sum(~inds),1) X(~inds,:)];
        yTest = Y(~inds);
%         xTest = standardize_test(xTest,sessInds(~inds),m,s);
        xTest =  standardize_test(xTest,ones(sum(~inds),1),m,s);
        
        
        % weight by percentage of positive and negative classes. I don't
        % totally get this but I am following Jim's code
        pos = mean(yTrain==1);
        neg = mean(yTrain==-1);
        mean_tmp = mean([pos neg]);
        pos = sprintf('%f',pos/mean_tmp);
        neg = sprintf('%f',neg/mean_tmp);
        param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
        
        % train on this data with this cv
        model = trainFun(double(yTrain),sparse(xTrain),param);
%         models = [models;model.w];
        
        % predict with model
        [pred, acc, dec] = testFun(double(yTest),sparse(xTest),model,'-b 0 -q');
%         dec_values(ind) = dec(:,1);
%         labels(ind) = yTest;
%         preds(ind) = pred;
        dec_values = [dec_values; dec(:,model.Label==1)];
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
[cBest,tBest,encBest,gBest]=ind2sub(size(auc_pen),ind);
cBest = Cs(cBest);
encBest = types(encBest);
if strcmpi(normType,'svm')        
    gBest = Gs(gBest);
else
    gBest = [];
end


function [yProbs,yPreds,yTest,A,err,C,G,tBest,encBest,Cs,Gs,Ts,aucs] = doRegFun(X,Y,T,folds,iFold,C,G,normType,sessInds,trialType,prctileThresh,nestedCvGroup,percentCV,useKfold,Cs,Gs,tBest,encBest,doCalcPenalty,doPermute)
% This does the classification.
% X = # trials x # features
% Y = # trials vector of responses
% folds = described above
% iFold = current fold number
% C = penalty parameter, if given
% normtype = L1 or L2
% sessInds = vector, same number of rows as X, identify the session for the
% row
%


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

    
% if no C given, calculate C for this fold.
if doCalcPenalty
    if ~useKfold
        [~,subFolds] = createFolds(sessInds(trainInds)',nestedCvGroup(trainInds));
    else
        [~,subFolds] = createKfolds(sum(trainInds),percentCV);
    end
    [C,G,tBest,encBest,Cs,Gs,Ts,aucs] = calcPenalty(xTrain,yTrainBool,T,subFolds,sessions,trialType,normType,prctileThresh,Cs,Gs);
else    
    Ts = tBest;
    aucs = NaN;
end
    
if encBest == 3
    trialType(:) = 3;
end

% set parameters for either L1 or L2
if strcmpi(normType,'L1')
    liblin_param = ['-c ' sprintf('%f',C) ' -s 6 -q'];
    trainFun = @train;
    testFun  = @predict;
    opts = '-b 1 -q';
elseif strcmpi(normType,'L2')
    liblin_param = ['-c ' sprintf('%f',C) ' -s 0 -q'];
    trainFun = @train;
    testFun  = @predict;
    opts = '-b 1 -q';
elseif strcmpi(normType,'svm')
    liblin_param = ['-t 2 -g ' sprintf('%f',G) ' -c ' sprintf('%f',C) ' -b 0 -q'];    
    trainFun = @svmtrain;
    testFun  = @svmpredict;
    opts = '-b 0 -q';
end


% make sure the classes are balanced by doing a number of subsampling
% rounds. First see how many observations to remove to balance them
numToRemove = sum(yTrainBool) - sum(~yTrainBool);
maxSamp = 1;
numToRemove = 0;
% numToRemove = 0;
% if numToRemove == 0
%     maxSamp = 1;
% end

% will store output from each round
preds = NaN(sum(~trainInds),maxSamp);
probs = NaN(sum(~trainInds),maxSamp);
% preds = NaN(1,maxSamp);
% probs = NaN(1,maxSamp);
 
for nSamp = 1:maxSamp
    
    % pick obs to remove
    toRemove = [];
    if numToRemove > 0
%         toRemove = randsample(find(yTrainBool==1),abs(numToRemove));
%         toRemove = randsample(find(Y(trainInds)==1),abs(numToRemove));
        toRemove = randsample(unique([find(~trainInds) find(Y==1)']),abs(numToRemove));
    elseif numToRemove < 0
%         toRemove = randsample(find(yTrainBool~=1),abs(numToRemove));
%         toRemove = randsample(find(Y(trainInds)~=1),abs(numToRemove));
        toRemove = randsample(unique([find(~trainInds) find(Y~=1)']),abs(numToRemove));
    end         
    trainIndsSamp = setdiff(find(trainInds),toRemove);
     
    % remove from training set
    % training set x
    xTrain = X(trainIndsSamp,T==tBest & trialType==encBest);
     
    % training set y
    yTrain = double(Y(trainIndsSamp));
    yTrain(yTrain==0) = -1;
    
     
    % standardize training data
%     [xTrain,m,s] = standardize(xTrain,sessInds(trainIndsSamp));
    [xTrain,m,s] = standardize(xTrain,ones(size(xTrain,1),1));
    
    % weight by class probabilities
    pos = mean(yTrain==1);
    neg = mean(yTrain==-1);
    mean_tmp = mean([pos neg]);
    pos = sprintf('%f',pos/mean_tmp);
    neg = sprintf('%f',neg/mean_tmp);
    param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
    
    % train model
    model = trainFun(double(yTrain),sparse(xTrain),param);
    if ~strcmpi(normType,'svm')
        if nSamp == 1        
            Ws = NaN(maxSamp,size(model.w,2));
        end
        Ws(nSamp,:) = model.w;
    else
        Ws = model.SVs;
    end
    
    % testing set
%     trainInds(toRemove) = 1;
    xTest = X(~trainInds,T==tBest & trialType==encBest);
    yTest = double(Y(~trainInds));
    yTest(yTest==0) = -1;
    
    % standardize testing data
%     xTest = standardize_test(xTest,sessInds(~trainInds),m,s);
    xTest = standardize_test(xTest,ones(sum(~trainInds),1),m,s);
    
    
%     probs = [];
%     fakeX = NaN(1,size(xTrain,2));
%     for i = 1:100
%     for c = 1:size(xTrain,2)
%        
%         r = randi(size(xTrain,1));
%         fakeX(c) = xTrain(r,c);
%         
%     end
%     [predSamp, acc, probSamp] = predict(double(yTest),sparse(fakeX),model,'-b 1 -q')
%     probs(i) = probSamp(1);
%     end
% %     % predict
% %     fakeX = randn(size(xTest));
% %     
% % %     fakeX = ones(size(xTest));
% % 	r = randi(length(trainIndsSamp))
% %     
%     ps = [];
%     accs = [];
%     prob = [];
%     for r = 1:size(xTrain,1)
%         fakeX = xTrain(r,:);
%         fakeY = yTrain(r);
% %         fakeX = X(trainIndsSamp(r),T==tBest & trialType==encBest);
% %         fakeX = standardize_test(fakeX,sessInds(~trainInds),m,s);        
%         [predSamp, acc, probSamp] = predict(double(yTest),sparse(fakeX),model,'-b 1 -q');
%         [predSamp, acc, probSamp] = predict(double(fakeY),sparse(fakeX),model,'-b 1 -q');
%         ps = [ps;predSamp];
%         accs = [accs;acc(1)];
%         prob = [prob;probSamp(1)];
% end
%     
%     [predSamp, acc, probSamp] = predict(double(yTest),sparse(xTest),model,'-b 1 -q');
%     
%     y = Y(toRemove);
%     y = double(y);
%     y(y==0) = -1;
%     x = X(toRemove,T==tBest & trialType==encBest);
%     x = standardize_test(x,sessInds(~trainInds),m,s); 
%      [predSamp, acc, probSamp] = predict(y,sparse(x),model,'-b 1 -q');
    % same as p = glmval(model.w',xTest,'logit','constant','off')?
    
    [predSamp, acc, probSamp] = testFun(double(yTest),sparse(xTest),model,opts);
    preds(:,nSamp) = predSamp;
    probs(:,nSamp) = probSamp(:,model.Label==1);
    
       
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



% tal = sessData.features.elecs(elecs);
function [powerData,sessInds,tal,powerParams] = loadAllPower(tal,subj,events,timeBins,eventsToUse,params)

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
powerData         = NaN(nEvents,nElecs,2,nTimes,nFreqBins);

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


% make encoding period one and two  a new dimension of matrix
encPeriods = [events.withinItemCount];
powEnc1    = powerDataAllFreqs(eventsToUse&encPeriods==1,:,:,:);
powEnc2    = powerDataAllFreqs(eventsToUse&encPeriods==2,:,:,:);
powerData  = cat(5,powEnc1,powEnc2);
powerData  = permute(powerData,[1,2,5,3,4]);


function [powerNorm,T,E,m,sigma,sessMap] = normPow_andMeanFreqs(powerData,sessions,params,freqs,m,sigma,sessMap)

[nEvents,nElecs,nEncs,nTimes,nFreqs] = size(powerData);
nFreqBins = size(params.freqBins,1);
powerNorm = NaN(nEvents,nElecs,nEncs,nTimes,nFreqBins);

% create matrix to keep track of time bins
timeMat = NaN(1,nElecs,nEncs,nTimes,nFreqBins);
for t = 1:size(timeMat,4)
    timeMat(:,:,:,t,:) = t;
end

% create matrix to keep track of enc periods
encMat = NaN(1,nElecs,nEncs,nTimes,nFreqBins);
for e = 1:size(encMat,3)
    encMat(:,:,e,:,:) = e;
end

% reshape into number of events by number of overall features for easier
% zscoring. I'm scoring here to put the different freqs on equal footing
powerDataTmp = reshape(powerData,size(powerData,1),[]);
uniqSess     = unique(sessions);
if isempty(m) || isempty(sigma)
    m       = NaN(length(uniqSess),size(powerDataTmp,2));
    sigma   = NaN(length(uniqSess),size(powerDataTmp,2));
    sessMap = NaN(length(uniqSess),1);
end

for s = 1:length(uniqSess)    
    sessInds = sessions==uniqSess(s);
    if all(isnan(m(s,:)))
        [powerDataTmp(sessInds,:),m(s,:),sigma(s,:)] = zscore(powerDataTmp(sessInds,:));
        sessMap(s) = uniqSess(s);
    else
        sessMapInd = sessMap == uniqSess(s);
        powerDataTmp(sessInds,:) = bsxfun(@rdivide, bsxfun(@minus, powerDataTmp(sessInds,:), m(sessMapInd,:)), sigma(sessMapInd,:));   
    end
end
    
% reshape back and bin by frequency
powerDataTmp = reshape(powerDataTmp,[nEvents,nElecs,nEncs,nTimes,nFreqs]);
if nFreqBins ~= 0    
    for f = 1:nFreqBins
        fInds = freqs >= params.freqBins(f,1)-.001 & freqs <= params.freqBins(f,2)+.001;
        powerNorm(:,:,:,:,f) = nanmean(powerDataTmp(:,:,:,:,fInds),5);
    end
else
    powerNorm = powerDataTmp;
end

% reshape into obs x features
powerNorm = reshape(powerNorm,nEvents,[]);
T         = reshape(timeMat,[1,nElecs*nEncs*nTimes*nFreqBins]);
E         = reshape(encMat,[1,nElecs*nEncs*nTimes*nFreqBins]);

% do I need to normalize again?? 


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














