function [auc,params] = TH1_refactor_phase(subj,params,saveDir)
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
auc = NaN;

% do we overwrite results file if it exists?
fname    = fullfile(saveDir,[subj '_class_pow.mat']);
powDir = fullfile(saveDir,'power');        
if ~exist(powDir,'dir');
    mkdir(powDir);
end

% if using power features
powFile  = fullfile(powDir,[subj '_power.mat']);
dataFunc = @loadAllPower;

% or if using phase features
if params.usePhase==1
    fname    = fullfile(saveDir,[subj '_class_phase.mat']);
    dataFunc = @loadAllPhase;
    powFile  = fullfile(powDir,[subj '_phase.mat']);
    if size(params.timeBins,1) > 1 || size(params.timeBins,2) > 1    
        fprintf('Please enter a single time point to use for phase.\n')
        return
    end
    
% or using both
elseif params.usePhase==2
    fname    = fullfile(saveDir,[subj '_class_phaseAndPow.mat']);
    dataFunc = @loadPowAndPhase;
    powFile  = fullfile(powDir,[subj '_phaseAndPow.mat']);
    if size(params.timeBins,1) > 1 && (params.timeBins(1) ~= params.timeBins(2))
        fprintf('Please enter a single time point to use for phase and power.\n')
        return
    end
    
% or phase lag!?    
elseif params.usePhase==3
    fname    = fullfile(saveDir,[subj '_class_phaseLag.mat']);
    dataFunc = @loadPhaseLag;
    powFile  = fullfile(powDir,[subj '_phaseLag.mat']); 
    
% power and phase lag    
elseif params.usePhase==4
    fname    = fullfile(saveDir,[subj '_class_powAndphaseLag.mat']);
    dataFunc = @loadPowerAndPhaseLag;
    powFile  = fullfile(powDir,[subj '_powAndPhaseLag.mat']); 
end

if exist(fname,'file') && params.saveOutput && ~params.overwrite
    fprintf('Classifier file already exists for %s.\n',subj)
    load(fname)
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
    useKfold      = params.useKfold;
    cvField       = params.cvField;
    nestedCvField = params.nestedCvField;
    k             = params.k;
    
     % load events
    if isfield(params,'diffRegress') && params.diffRegress == 1
        events = load(fullfile('/scratch/jfm2/diffRegress/RAM_TH1/',[subj '_events.mat']));
        events = events.events;
    else
        events = get_sub_events('RAM_TH1',subj);                
    end
    events = addTranspose_distErr(events);
    
    % filter to just specific confidence levels
    if ~isfield(params,'confsToUse')
        confsToUse = 0:2;
    else
        confsToUse = params.confsToUse;
    end    
    eventsToUse = params.eventFilter(events) & ismember([events.confidence],confsToUse);   
    
    % limit correct items to high conf, limit incorrect to not-high conf
    if isfield(params,'strongRecVsWeak')
        strongAndWeakRec = ([events.confidence]==2 & [events.recalled]==1) | ([events.confidence]<2 & [events.recalled]==0);        
        eventsToUse      = eventsToUse & strongAndWeakRec;        
    end
    
    
    % if not enough events in a session, exclude
    sessions  = [events.session];
    evPerSess = grpstats(sessions(eventsToUse),sessions(eventsToUse),'numel');
    uniqSess  = unique(sessions);
    for s = 1:length(uniqSess)
        if evPerSess(s) < 15
            eventsToUse(sessions==uniqSess(s)) = 0;
        end
    end
    
    % new session vector
    sessions  = [events(eventsToUse).session];
    
    % if we have no events, or we don't have enough sessions for the
    % selector type, return
    doFolds = 1;
    if sum(eventsToUse) == 0
        fprintf('No enough events for %s.\n',subj)
        return
    elseif strcmp(cvField,'session')
        if length(unique(sessions)) == 2 && strcmp(nestedCvField,'session')
            fprintf('No enough sessions inner and outer session CV for %s. ',subj)
            fprintf('Creating full model only.\n')
            doFolds = 0;
        elseif length(unique(sessions)) == 1% && strcmp(nestedCvField,'session')
            params.cvField = 'trial';
            cvField = 'trial';
%             fprintf('No enough sessions inner and outer session CV or full model for %s.\n',subj)
            fprintf('No enough sessions session CV for %s. Switching to trial.\n',subj)
%             return
        end
    end
               
    % outer fold selector
    innerSel   = [events(eventsToUse).(nestedCvField)];
    outerSel   = [events(eventsToUse).(cvField)];
    
    % response data    
    distanceErrs = [events(eventsToUse).distErr]';
    confs        = [events(eventsToUse).confidence]';
    if isfield(params,'classifyConf') && params.classifyConf
        Y = confs == 2;
    else
        if isfield(params,'useMedianAsThresh') && params.useMedianAsThresh
            params.correctThresh = median(distanceErrs);
        end
        if params.correctThresh > 0
            Y = distanceErrs <= params.correctThresh;
        else
            if isfield(params,'countTransposedAsCorrect') && params.countTransposedAsCorrect
                Y = [[events(eventsToUse).recalled]==1 | [events(eventsToUse).recalled_ifFlipped]==1]';
            else
                Y = [events(eventsToUse).recalled]';
            end            
        end
        if isfield(params,'countLowConfAsWrong') && params.countLowConfAsWrong
           Y(confs==0) = 0; 
        end
    end
    
    if doPermute
        randOrder = randperm(length(Y));
        Y = Y(randOrder);
        distanceErrs = distanceErrs(randOrder);
        confs = confs(randOrder);
    end
    
    % load power
    if ~params.loadPower
        if ~isfield(params,'doBipol')
            tal = getBipolarSubjElecs(subj,1,1,params.excludeEpiElecs);
        else
            tal = getBipolarSubjElecs(subj,params.doBipol,1,params.excludeEpiElecs);
        end
        
        % filter out bad chans using Andrew's method
        if isfield(params,'excludeBadChansWatrous') && params.excludeBadChansWatrous
            uniqSess     = unique(sessions);
            badChans     = [];
            for s = 1:length(uniqSess)                
                badChanFile = sprintf('/home2/andrew.watrous/Preproc/RAM_TH1/Screened/%s_sess%d.mat',subj,uniqSess(s));
                if ~exist(badChanFile,'file')
                    fprintf('excludeBadChansWatrous requested but bad channel file %s not found. Skipping %s.\n',badChanFile,subj)
                    return
                end
                badChanInfo = load(badChanFile);
                badChanInds = badChanInfo.screened_channels.bad_electrodes;
                badChans    = [badChans;vertcat(badChanInfo.screened_channels.elec_info(badChanInds==1).channel)];
            end
            badChans = unique(badChans,'rows');
            fprintf('Excluding %d channels from %s.\n',size(badChans,1),subj)
            
            % remove from tal
            tal = tal(~ismember(vertcat(tal.channel),badChans,'rows'));            
        end
        [powMat,timeMat,~,tal,powParams] = dataFunc(tal,subj,events,params.timeBins,eventsToUse,params);        
        if params.savePower
            save(powFile,'powMat','timeMat','tal','powParams','sessions')
        end
    else
        load(powFile);
    end        
    
    % create outer cross val folds
    if useKfold
        folds = create_stratKfolds(Y,k);
        nFolds = k;
    else
        [nFolds,folds,trials] = createFolds(sessions,outerSel);
    end
    
    res = [];
    [res.yProbs, res.yPreds, res.yTests, res.bestTime, res.bestC, res.aucFold] = deal(cell(1,nFolds));
    if doFolds
        for thisFold = 1:nFolds
            fprintf('%s: fold %d of %d.\n',subj,thisFold,nFolds)
            [yProbFold,yPredFold,yTestFold,tFold,cFold] = run_fold_with_nested_folds(powMat,timeMat,Y,sessions,folds,thisFold,innerSel,params,powParams);
            
            % store outputs from fold
            res.yProbs{thisFold} = yProbFold;
            res.yPreds{thisFold} = yPredFold;
            res.yTests{thisFold} = yTestFold;            
            res.bestTime{thisFold} = tFold;
            res.bestC{thisFold}    = cFold;
            
            if ismember([-1,1],vertcat(res.yTests{:}))
                [~,~,~,aucTmp] = perfcurve(vertcat(res.yTests{:}),vertcat(res.yProbs{:}),1);
                fprintf('%s: running AUC = %.3f\n',subj,aucTmp)
            end
            
            if params.useKfold
                [~,~,~,res.aucFold{thisFold}] = perfcurve(yTestFold,yProbFold,1);
                fprintf('%s: running mean fold AUC = %.3f\n',subj,mean(vertcat(res.aucFold{:})))
            end
        end
        auc = aucTmp;
        res.auc = aucTmp;
    end
    
    % make full model
    fprintf('%s: Creating full model.\n',subj)
    res.fullModel = [];
    [model,thisT,thisC] = make_full_model(powMat,timeMat,Y,sessions,folds,params,powParams);
    res.fullModel.model = model;    
    res.fullModel.bestT = thisT;
    res.fullModel.bestC = thisC;    
    
    if saveOutput
        save(fname,'res','params','tal','distanceErrs','eventsToUse','cvField','confs')
    end
    

catch e
    fname = fullfile(saveDir,[subj '_class_error.mat']);
    save(fname,'e')
end

function [model,thisT,thisC] = make_full_model(X,T,Y,sessions,folds,params,powParams)
% convert Y into double for liblinear
Y = double(Y);
Y(Y==0) = -1;

% unique times and enc periods
Ts = 1:size(params.timeBins,1);
Cs = params.Cs;

% run nested fold to find optimal time bin, and penalty
if length(Ts) ~= 1 || length(Cs) ~= 1
    [aucs] = run_nested_fold(X,T,Y,sessions,folds,params,powParams);
    ind=find(aucs==max(aucs(:)),1,'first');
    [tBestInd,cBestInd]=ind2sub(size(aucs),ind);
    thisT = Ts(tBestInd);
    thisC = Cs(cBestInd);
    fprintf('Max AUC = %.3f, best time = %d, best C = %.4f\n',max(aucs(:)),thisT,thisC)
else
    thisT = Ts;
    thisC = Cs;
end
% balance the classes
% X(toRemove,:) = [];
% Y(toRemove) = [];
% sessions(toRemove) = [];

% normalize full training set        
[X,m,s,sessMap] = normPow(X,sessions,[],[],[]);

% mask encoding period and time bin
featureMask = T == thisT;
X = X(:,featureMask);

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
pos = mean(Y==1);
neg = mean(Y==-1);
mean_tmp = mean([pos neg]);
pos = sprintf('%f',pos/mean_tmp);
neg = sprintf('%f',neg/mean_tmp);
param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
% param = liblin_param;
model = trainFun(double(Y),sparse(X),param);


function [yProb,yPred,yTest,thisT,thisC] = run_fold_with_nested_folds(X,T,Y,sessions,folds,thisFold,innerSel,params,powParams)

% convert Y into double for liblinear
Y = double(Y);
Y(Y==0) = -1;

% unique times and enc periods
Ts = 1:size(params.timeBins,1);
Cs = params.Cs;

% get training data
trainMask = folds(thisFold,:);
if params.useKfold
    subFolds = create_stratKfolds(Y(trainMask),params.k);    
else
    [~,subFolds] = createFolds(sessions(trainMask),innerSel(trainMask));
end
xTrain    = X(trainMask,:);
yTrain    = Y(trainMask);
sessTrain = sessions(trainMask);

% run nested fold to find optimal enc, time bin, and penalty
if length(Ts) ~= 1 || length(Cs) ~= 1
    [aucs] = run_nested_fold(xTrain,T,yTrain,sessTrain,subFolds,params,powParams);
    ind=find(aucs==max(aucs(:)),1,'first');
    [tBestInd,cBestInd]=ind2sub(size(aucs),ind);
    thisT = Ts(tBestInd);
    thisC = Cs(cBestInd);
    fprintf('Max AUC = %.3f, best time = %d, best C = %.4f\n',max(aucs(:)),thisT,thisC)
else
    thisT = Ts;
    thisC = Cs;
end

% balance the classes if not doing statified k-fold cv
% if ~params.useKfold && strcmp(params.nestedCvField,'blocknum')    
%     pos = find(yTrain == 1);
%     neg = find(yTrain == -1);
%     
%     numToRemove = length(pos)-length(neg);
%     toRemove = [];
%     if numToRemove > 0
%         toRemove=randsample(pos,abs(numToRemove));
%     elseif numToRemove < 0
%         toRemove=randsample(neg,abs(numToRemove));
%     end
%     
%     % balance the classes
%     xTrain(toRemove,:) = [];
%     yTrain(toRemove) = [];
%     sessTrain(toRemove) = [];
% end

% normalize full training set        
[xTrain,m,s,sessMap] = normPow(xTrain,sessTrain,[],[],[]);
% [xTrain,T,E,m,s,sessMap] = normPow_andMeanFreqs(xTrain,sessTrain,params,powParams.params.pow.freqs,[],[],[]);

% mask encoding period and time bin
featureMask = T == thisT;
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
pos = mean(yTrain==1);
neg = mean(yTrain==-1);
mean_tmp = mean([pos neg]);
pos = sprintf('%f',pos/mean_tmp);
neg = sprintf('%f',neg/mean_tmp);
param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
% param = liblin_param;
model = trainFun(double(yTrain),sparse(xTrain),param);

% test data
xTest    = X(~trainMask,:);
yTest    = Y(~trainMask);
sessTest = sessions(~trainMask);

% normalize test data and mean into freqs of interest
if strcmp(params.cvField,'session')
    xTest = normPow(xTest,sessTest,[],[],[]);
else
    xTest = normPow(xTest,sessTest,m,s,sessMap);
end
% xTest = normPow_andMeanFreqs(xTest,sessTest,params,powParams.params.pow.freqs,m,s,sessMap);

% test on held out
xTest = xTest(:,featureMask);
[pred, acc, dec] = testFun(double(yTest),sparse(xTest),model,'-b 1 -q');

% probability of positive class
yProb = dec(:,model.Label==1);
yPred = pred;

function [aucs] = run_nested_fold(X,T,Y,sessions,subFolds,params,powParams)

% convert Y into double for liblinear
Y = double(Y);
Y(Y==0) = -1;

% balance the classes if not doing statified k-fold cv
subFolds = double(subFolds);
if ~params.useKfold && strcmp(params.nestedCvField,'blocknum')    
    for thisFold = 1:size(subFolds,1)
        pos = find(subFolds(thisFold,:)==1 & Y' == 1);
        neg = find(subFolds(thisFold,:)==1 & Y' == -1);

        numToRemove = length(pos)-length(neg);    
        toRemove = [];
        if numToRemove > 0
            toRemove=randsample(pos,abs(numToRemove));
        elseif numToRemove < 0
            toRemove=randsample(neg,abs(numToRemove));
        end    
        subFolds(thisFold,toRemove) = NaN;
    end
end

% unique times and enc periods
Ts = 1:size(params.timeBins,1);
Cs = params.Cs;

% will hold auc for each combination of parametres
aucs = NaN(length(Ts),length(Cs));
probs = NaN(length(Ts),length(Cs),size(subFolds,2));

% max loop number for parfor
maxInd = length(Cs)*length(Ts);

% parallelize over which ever is greater: number of folds or number or
% elements to grid search over
if size(subFolds,1) > maxInd && strcmp(params.nestedCvField,'blocknum') && ~params.useKfold
    for thisFold = 1:size(subFolds,1)
        
        trainMask  = subFolds(thisFold,:)==1;
        testMask   = subFolds(thisFold,:)==0;
        
        % training data for fold
        xTrain    = X(trainMask,:);
        yTrain    = Y(trainMask);
        sessTrain = sessions(trainMask);
        
        % normalize training data and mean into frequencies of interest
        [xTrain,m,s,sessMap] = normPow(xTrain,sessTrain,[],[],[]);
        
        % test data
        xTest    = X(testMask,:);
        yTest    = Y(testMask);
        sessTest = sessions(testMask);
        
        % normalize test data and mean into freqs of interest
        if strcmp(params.nestedCvField,'session')
            xTest = normPow(xTest,sessTest,[],[],[]);
        else
            xTest = normPow(xTest,sessTest,m,s,sessMap);
        end
                
        probs_fold = NaN(length(Ts),length(Cs),sum(testMask));
        for ind = 1:maxInd
            
            % get current, E, T, C
            [t,c] = ind2sub([length(Ts),length(Cs)],ind);            
            thisT = Ts(t);
            thisC = Cs(c);
            
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
            featureMask = T == thisT;
            xTrain_fold = xTrain(:,featureMask);
            
            % create model
            param = liblin_param;
            model = trainFun(double(yTrain),sparse(xTrain_fold),param);
            
            % test on held out
            xTest_fold = xTest(:,featureMask);
            [pred, acc, dec] = testFun(double(yTest),sparse(xTest_fold),model,'-b 1 -q');
            probs_fold(t,c,:) = dec(:,1);
        end
        
        % this will only work with leave one item out (blocknum cv) because
        % parfor loops are the fucking worst. i can't index it correctly
        probs(:,:,thisFold) = probs_fold;
        
    end
    
    % compute auc for each encoding period, time bin, C parameter    
    for t = 1:length(Ts)
        for c = 1:length(Cs)
            [~,~,~,aucs(t,c)] = perfcurve(Y,squeeze(probs(t,c,:)),1);
        end
    end
    
    
% if we are parallel over hyper parameters instead of folds
else
    for ind = 1:maxInd
        
        % get current, E, T, C
        [thisT,thisC] = ind2sub([length(Ts),length(Cs)],ind);
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
            
            % masks 
            trainMask  = subFolds(thisFold,:)==1;
            testMask   = subFolds(thisFold,:)==0;
            
            % training data for fold
            xTrain    = X(trainMask,:);
            yTrain    = Y(trainMask);
            sessTrain = sessions(trainMask);
            
            % normalize training data and mean into frequencies of interest
            [xTrain,m,s,sessMap] = normPow(xTrain,sessTrain,[],[],[]);
            %         [xTrain,T,E,m,s,sessMap] = normPow_andMeanFreqs(xTrain,sessTrain,params,powParams.params.pow.freqs,[],[],[]);
            
            % mask encoding period and time bin
            featureMask = T == thisT;
            xTrain = xTrain(:,featureMask);
            
            % weight by percentage of positive and negative classes
            %         pos = mean(yTrain==1);
            %         neg = mean(yTrain==-1);
            %         mean_tmp = mean([pos neg]);
            %         pos = sprintf('%f',pos/mean_tmp);
            %         neg = sprintf('%f',neg/mean_tmp);
            %         param = [liblin_param ' -w1 ' num2str(pos) ' -w-1 ' num2str(neg)];
            
            % create model
            param = liblin_param;
            model = trainFun(double(yTrain),sparse(xTrain),param);
            
            % test data
            xTest    = X(testMask,:);
            yTest    = Y(testMask);
            sessTest = sessions(testMask);
            
            % normalize test data and mean into freqs of interest
            if strcmp(params.nestedCvField,'session')
                xTest = normPow(xTest,sessTest,[],[],[]);
            else
                xTest = normPow(xTest,sessTest,m,s,sessMap);
            end
            %         xTest = normPow_andMeanFreqs(xTest,sessTest,params,powParams.params.pow.freqs,m,s,sessMap);
            
            % test on held out
            xTest = xTest(:,featureMask);
            [pred, acc, dec] = testFun(double(yTest),sparse(xTest),model,'-b 1 -q');
            
            dec_values = [dec_values; dec(:,model.Label==1)];
            labels = [labels;yTest];
        end
        [~,~,~,aucs(ind)] = perfcurve(labels,dec_values,1);
    end
end


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
function [powerData,timeMat,sessInds,tal,powerParams] = loadAllPower(tal,subj,events,timeBins,eventsToUse,params)

% power matrix will initially be # events x # elecs x # times x # freqs
nTimes  = size(timeBins,1);
nEvents = sum(eventsToUse);
[tal]   = filterTalByRegion(tal,params.region);
nElecs  = length(tal);

% load parameters used to create power
powerParams = load(fullfile(params.powerPath,'params_RAM_TH1.mat'));
freqs       = powerParams.params.pow.freqs;
nFreqs      = length(powerParams.params.pow.freqs);
nFreqBins   = size(params.freqBins,1);

% will hold power data
powerDataAllFreqs = NaN(length(events),nElecs,nTimes,nFreqs);
powerData         = NaN(sum(eventsToUse),nElecs,nTimes,nFreqBins);

% 
basePath  = params.powerPath;
subjPath  = fullfile(basePath,subj);
sessions  = [events.session];
uniqSess  = unique(sessions);

% loop over each electrode
fprintf('%s: Loading power data for %d electrodes.\n',subj,nElecs);
for e = 1:nElecs
%     fprintf('%s: Loading elec %d of %d.\n',subj,e,nElecs);
    elecNum = tal(e).channel;
    
    % and session
    for s = 1:length(uniqSess)
        sessInds = sessions==uniqSess(s);        
        
        % load power for electrode/session
        if ~isfield(params,'doBipol') || params.doBipol
            fname = fullfile(subjPath,'RAM_TH1_events',num2str(uniqSess(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
        else
            fname = fullfile(subjPath,'RAM_TH1_events',num2str(uniqSess(s)),[num2str(elecNum(1)),'.mat']);
        end
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
        powerData(eventsToUse,:,:,f) = nanmean(powerDataAllFreqs(eventsToUse,:,:,fInds),4);
    end
else
    powerData = powerDataAllFreqs(eventsToUse,:,:,:);
end

% create matrix to keep track of time bins
if nFreqBins == 0;nFreqBins=nFreqs;end
timeMat = NaN(1,nElecs,nTimes,nFreqBins);
for t = 1:size(timeMat,3)
    timeMat(:,:,t,:) = t;
end
 
% reshape into obs x features
powerData = reshape(powerData,size(powerData,1),[]);
timeMat   = reshape(timeMat,[1,nElecs*nTimes*nFreqBins]);

function [phaseData,timeMat,sessInds,tal,powerParams] = loadAllPhase(tal,subj,events,timeBin,eventsToUse,params,returnRaw)

if ~exist('returnRaw','var') || isempty(returnRaw)
    returnRaw = 0;
end

% power matrix will initially be # events x # elecs x # times x # freqs
nTimes  = length(timeBin(1):timeBin(end));
nEvents = sum(eventsToUse);
[tal]   = filterTalByRegion(tal,params.region);
nElecs  = length(tal);

% load parameters used to create power
powerParams = load(fullfile(params.powerPath,'params_RAM_TH1.mat'));
nFreqs      = length(powerParams.params.pow.freqs);

% will hold phase data
phaseData = NaN(length(events),nElecs,nTimes,nFreqs);


basePath  = params.powerPath;
subjPath  = fullfile(basePath,subj);
sessions  = [events.session];
uniqSess  = unique(sessions);

% loop over each electrode
fprintf('%s: Loading phase data for %d electrodes.\n',subj,nElecs);
for e = 1:nElecs
%     fprintf('%s: Loading elec %d of %d.\n',subj,e,nElecs);
    elecNum = tal(e).channel;
    
    % and session
    for s = 1:length(uniqSess)
        sessInds = sessions==uniqSess(s);        
        
        % load power for electrode/session
        if ~isfield(params,'doBipol') || params.doBipol
            fname = fullfile(subjPath,'RAM_TH1_events',num2str(uniqSess(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
        else
            fname = fullfile(subjPath,'RAM_TH1_events',num2str(uniqSess(s)),[num2str(elecNum(1)),'.mat']);
        end
        sessPhase = load(fname);
        
        elecPhase = sessPhase.sessOutput.phase(:,timeBin,:);      
        elecPhase = permute(elecPhase,[3 2 1]);
        phaseData(sessInds,e,:,:) = elecPhase;
        
    end
end

phaseData    = phaseData(eventsToUse,:,:,:);
nPhaseObs    = 1;
if ~returnRaw
    phaseDataSin = sin(phaseData);
    phaseDataCos = cos(phaseData);
    phaseData    = cat(5,phaseDataSin,phaseDataCos);
    nPhaseObs    = 2;
end

% create matrix to keep track of time bins
timeMat = NaN(1,nElecs,nTimes,nFreqs,nPhaseObs);
for t = 1:size(timeMat,3)
    timeMat(:,:,t,:,:) = t;
end
 
% reshape into obs x features
if ~returnRaw
    phaseData = reshape(phaseData,size(phaseData,1),[]);
    timeMat   = reshape(timeMat,[1,nElecs*nTimes*nFreqs*nPhaseObs]);
end
function [powPhaseData,timeMat,sessInds,tal,powerParams] = loadPowAndPhase(tal,subj,events,timeBin,eventsToUse,params)
% this is kind of ineffecient because I'm loading each file twice

[powerData,timeMat,sessInds,tal,powerParams] = loadAllPower(tal,subj,events,timeBin,eventsToUse,params);
[phaseData,timeMatPhase] = loadAllPhase(tal,subj,events,timeBin(1),eventsToUse,params); 
powPhaseData = cat(2,powerData,phaseData);
timeMat      = cat(2,timeMat,timeMatPhase);

function [powPhaseLagData,timeMat,sessInds,tal,powerParams] = loadPowerAndPhaseLag(tal,subj,events,timeBin,eventsToUse,params)
% % this is kind of ineffecient because I'm loading each file twice
% 
[powerData,timeMat,sessInds,tal,powerParams] = loadAllPower(tal,subj,events,timeBin,eventsToUse,params);
[phaseLagData,timeMatPhaseLag] = loadPhaseLag(tal,subj,events,timeBin(1),eventsToUse,params);
powPhaseLagData = cat(2,powerData,phaseLagData);
timeMat      = cat(2,timeMat,timeMatPhaseLag);

function [phaseLag,timeMat,sessInds,tal,powerParams] = loadPhaseLag(tal,subj,events,timeBin,eventsToUse,params)

[phaseData,~,sessInds,tal,powerParams] = loadAllPhase(tal,subj,events,timeBin(1):timeBin(end),eventsToUse,params,1); 

% compute all pairs of phase lags
nElecs   = size(phaseData,2);
allPairs = nchoosek(1:nElecs,2);
fprintf('%s: Computing pairwise phase lags for %d pairs.\n',subj,length(allPairs));
phaseLag = NaN(size(phaseData,1),length(allPairs),1,size(phaseData,4),'single');
for pair = 1:length(allPairs)
    phaseLagPair = phaseData(:,allPairs(pair,1),:,:)-phaseData(:,allPairs(pair,2),:,:);
    phaseLagPair = mod(phaseLagPair+pi,2*pi)-pi;    
    phaseLag(:,pair,:,:) = circ_mean(phaseLagPair,[],3);
%     phaseLag(:,pair,:,:) = abs(mean(exp(1i*phaseLagPair),3));
end

% mean by time
% phaseLag = nanmean(phaseLag,3);

timeMat = NaN(1,length(allPairs),size(phaseLag,3),size(phaseLag,4));
for t = 1:size(timeMat,3)
    timeMat(:,:,t,:) = t;
end

% reshape into obs x features
phaseLag = reshape(phaseLag,size(phaseLag,1),[]);
timeMat   = reshape(timeMat,1,[]);

% powPhaseData = cat(2,powerData,phaseData);
% timeMat      = cat(2,timeMat,timeMatPhase);

%mod(a-b+pi,2*pi)-pi
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
    

function events = addTranspose_distErr(events)
xCenter = 384.8549;
yCenter = 358.834;
x = [events.chosenLocationX]-xCenter;
y = [events.chosenLocationY]-yCenter;
xChest = [events.locationX]-xCenter;
yChest = [events.locationY]-yCenter;

distErr_transpose = sqrt((-xChest - x).^2 + (-yChest - y).^2);
tmp=num2cell(distErr_transpose);
[events.distErr_transpose] = deal(tmp{:});
recalled_ifFlipped = NaN(1,length(events));
recalled_ifFlipped([events.isRecFromStartSide]==0) = distErr_transpose([events.isRecFromStartSide]==0)<=13;
tmp=num2cell(recalled_ifFlipped);
[events.recalled_ifFlipped] = deal(tmp{:});












