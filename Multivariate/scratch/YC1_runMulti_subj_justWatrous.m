function [perf,AUC,r,subject,params,mse] = YC1_runMulti_subj(subj,params,saveDir)
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
mse     = [];
sse     = [];
mae     = [];
r2      = [];
r       = [];

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
    powParams = load(fullfile(params.powerPath,'params.mat'));
    
    % Setting time bins for convenience:
    tEnds     = (powParams.params.pow.timeWin:powParams.params.pow.timeStep:powParams.params.eeg.durationMS)+powParams.params.eeg.offsetMS;
    tStarts   = tEnds - powParams.params.pow.timeWin+1;
    powParams.timeBins = [tStarts' tEnds'];
    
    % load events
    events = get_sub_events('RAM_YC1',subj);
    
    % add the test error to the learning trials
    events  = addErrorField(events);
    session = [events.session];
    
    % filter to events of interest
    eventsToUse = params.eventFilter(events);
    if sum(eventsToUse) < 10
        fprintf('Not enough events for %s.\n',subj)
        return
    end
    
    % and filter to encoding period of interest
    [~,firstIdx,~] = unique([session(eventsToUse)' [events(eventsToUse).blocknum]'],'rows','first');
    if any(strcmp({'first','second'},params.encPeriod))        
        inds = find(eventsToUse);
        if strcmp(params.encPeriod,'first')
            eventsToUse(inds(firstIdx+1)) = false;
        elseif strcmp(params.encPeriod,'second')
            eventsToUse(inds(firstIdx+1)) = false;
        end
    end    
    
    % get parameters
    freqBins      = params.freqBins;
    timeBins      = params.timeBins;
    modelEachTime = params.modelEachTime;
    doBinary      = params.doBinary;
    saveOutput    = params.saveOutput;
    doPermute     = params.doPermute;
    normType      = params.normType;
    
    % load power for all electrodes
    powDir = fullfile(saveDir,'power');
    powFile = fullfile(powDir,[subj '_binnedPower.mat']);
    if ~exist(powDir,'dir')
        mkdir(powDir)
    end
    if params.loadPower
        powerData = load(powFile);
        powerData = powerData.powerData;
    else
        powerData = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params);
        powerData = permute(powerData,[3 1 2 4]);
        
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
        powerData = wstPPC;
                
        if params.savePower
            powFile = fullfile(powDir,[subj '_binnedPower.mat']);
            save(powFile,'powerData','params')
        end
    end
         
    
    % if params.encPeriod is 'combined', then reshape the matrix to add an
    % additional dimensions that is encoding trial number (1 or 2). So
    % feature matrix will be:
    %       nItems x nFreqs x nTimes x nElecs x 2    
    
    if strcmpi(params.encPeriod,'combined') || strcmpi(params.encPeriod,'average')        
        
        % this is assuming every other entry is encoding 1, encoding 2,
        % encoding 1, encoding 2. This *should* always be true
        firstEnc  = powerData(1:2:end,:,:,:);
        secondEnc = powerData(2:2:end,:,:,:);        
        powerData = cat(5,firstEnc,secondEnc);
        
        if strcmpi(params.encPeriod,'average')
            powerData = nanmean(powerData,5);
        end
        
    end
    
    
    
    % size of feature matrix
    nElecs = size(powerData,4);
    nTimes = size(powerData,3);
    nFreqs = size(powerData,2);    
    nItemObs = size(powerData,5);        
    
    % response data
    Y = [events(eventsToUse).testError]';
    if doBinary
        Y  = Y < median(Y);
    end    
    
    % determine the cross validation folds. each fold actually leaves out one
    % learning pair, not one trial. For example, if you had 5 learning pairs,
    % the folds would look like, where each row represents the hold outs and
    % the training for that fold:
    %
    % 0 0 1 1 1 1 1 1 1 1
    % 1 1 0 0 1 1 1 1 1 1
    % 1 1 1 1 0 0 1 1 1 1
    % 1 1 1 1 1 1 0 0 1 1
    % 1 1 1 1 1 1 1 1 0 0
    %
    % JFM edit 9/16/2015: if params.encPeriod is 'first' or 'second' or
    % 'combined', it will not hold out pairs because there are no pairs.
    % If params.encPeriod is 'both', it will hold out pairs. The below
    % logic is unchanged.
    %
    % Note: this is not influenced by params.nCV. params.nCV is the number
    % of cross validation folds to estimate lambda with lassoglm. For the
    % training test iterations, we are always doing leave one object out.    

    [trials,~,trialInds] = unique([session(eventsToUse)' [events(eventsToUse).blocknum]'],'rows');
    if strcmpi(params.encPeriod,'combined') || strcmpi(params.encPeriod,'average')
        trialInds = trialInds(firstIdx);    
        Y = Y(firstIdx);
    end
    
    nFolds = size(trials,1);
    folds = false(nFolds,size(trialInds,1));
    for iFold = 1:nFolds
        folds(iFold,:) = trialInds ~= iFold;
    end  
    
    % use hagai's error metric instead of the normal one?
    if isfield(params,'useHagai') && params.useHagai == 1
        clusteredY = load('/home1/jfm2/matlab/YC_analyses/Multivariate/clusteredY2.mat');
        sInd = strcmp(clusteredY.subjects,subj);
        hagaiLabel = clusteredY.binaryBehaviorMetric{sInd};
        Y_tmp = [hagaiLabel';hagaiLabel'];
        Y_tmp = Y_tmp(:);
        if length(Y) ~= length(Y_tmp)
            fprintf('problem with hagai labels.\n')
            return
        end
        Y = logical(Y_tmp-1);
    end
    
    % permute the responses if desired
    if doPermute
        randOrder = randperm(size(trials,1));        
        if ~isfield(params,'encPeriod') || strcmpi(params.encPeriod,'both')
            randOrder = [randOrder;randOrder];
            randOrder = randOrder(:);
            Ytmp = reshape(Y,2,[]);
            Y = reshape(Ytmp(:,randOrder),[],1);
        else
            Y = Y(randOrder);
        end        
    end
    

    objLocs = vertcat(events(eventsToUse).objLocs);
    
    % We can model time points seperately, so # features = # freqs x # elecs,
    % or we can model it all together, so # features = # times x # freqs x #
    % elecs.
    res = [];
    if modelEachTime
        perf = NaN(1,nTimes);
        lambdas = NaN(1,nTimes);
        for t = 1:nTimes
            
            % reshape into # trials x # features
            X = reshape(squeeze(powerData(:,:,t,:,:)),size(powerData,1),nFreqs*nElecs*nItemObs);
            
            % see if we are precomputing lambda based on all the data
            lambda = [];
            if params.crossValStrictness == 0
                if ~isempty(params.lambda)
                    lambda = params.lambda(t);
                else
                    if params.alpha < 1 && strcmp(normType,'L1')
                        fprintf('Subject %s: Computing optimal lambda using elastic net normalization.\n',subj)
                    else
                        fprintf('Subject %s: Computing optimal lambda using %s normalization.\n',subj,normType)
                    end
                    if strcmp(normType,'L1')
                        [stats,lambda] = calcLambda(X,Y,doBinary,params.nCV,params.alpha);
                    elseif strcmp(normType,'L2')
                        lambda = calcPenalty(X,Y,params.nCV);
                    end
                end
                lambdas(t) = lambda;
            end
            
            % will hold results from each fold
            [res(t).yPred,res(t).yTest,res(t).A,res(t).intercept,res(t).err] = deal(cell(nFolds,1));
            
            % run for each fold
            for iFold = 1:nFolds
                fprintf('Subject %s: Time %d of %d, Fold %d of %d.\n',subj,t,nTimes,iFold,nFolds)
                [res(t).yPred{iFold},...
                    res(t).yTest{iFold},...
                    res(t).A{iFold},...
                    res(t).intercept{iFold},...
                    res(t).err{iFold}] = doRegFun(X,Y,folds(iFold,:),lambda,params.nCV,params.alpha,normType);
            end
            
            
            if doBinary
                perf(t) = mean(vertcat(res(t).err{:}));
                res(t).perf = perf(t);
                yPred = vertcat(res(t).yPred{:});
                [~,~,~,res(t).AUC] = perfcurve(Y,yPred,true);
                %                res(t).AUC = compute_auc(yPred,Y);
                AUC(t) = res(t).AUC;
                
            else
                mse(t) = mean(vertcat(res(t).err{:}).^2);
                res(t).mse = mse(t);
                sse(t) = sum(vertcat(res(t).err{:}).^2);
                res(t).sse = sse(t);
                mae(t) = mean(abs(vertcat(res(t).err{:})));
                res(t).mae = mae(t);
                
                ssr = sum((vertcat(res(t).yTest{:})-vertcat(res(t).yPred{:})).^2);
                %                 sst = var(vertcat(res(t).yTest{:}))*length(vertcat(res(t).yTest{:}));
                sst = sum((vertcat(res(t).yTest{:}) - mean(vertcat(res(t).yTest{:}))).^2);
                
                r2(t)  = 1 - ssr/sst;
                res(t).r2 = r2(t);
                
                r(t) = corr(vertcat(res(t).yPred{:}),vertcat(res(t).yTest{:}));
                res(t).r = r(t);
                
                %      optimal_weights = w(:,FitInfo.IndexMinMSE);
                %  SStot = var(y)*length(y);
                %  predicted_values = X*optimal_weights;
                %  SSres = sum( (y(:)-predicted_values(:)).^2 );
                %  R2 = 1 - SSres/SStot;
                
            end
        end
        lambda = lambdas;
        
        % if using all time points in one model, current what it is set to do
    else
        
        % reshape into # trials x # features
        X = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs*nItemObs);
        
        % see if we are precomputing lambda based on all the data
        lambda = [];
        if params.crossValStrictness == 0
            if ~isempty(params.lambda)
                lambda = params.lambda;
            else
                
                if params.alpha < 1 && strcmp(normType,'L1')
                    fprintf(['Subject %s: Computing optimal lambda using elastic net normalization.\n'],subj)
                else
                    fprintf(['Subject %s: Computing optimal lambda using %s normalization.\n'],subj,normType)
                end
                if strcmp(normType,'L1')
                    [stats,lambda] = calcLambda(X,Y,doBinary, ...
                        params.nCV,params.alpha);
                elseif strcmp(normType,'L2')
                    lambda = calcPenalty(X,Y,params.nCV);
                end
                
                
                %                fprintf('Subject %s: Computing optimal lambda.\n',subj)
                %                [stats,lambda] = calcLambda(X,Y,doBinary,params.nCV,params.alpha);
            end
        end
        
        % run for each fold
        [res.yPred,res.yTest,res.A,res.intercept,res.err] = deal(cell(nFolds,1));
        for iFold = 1:nFolds
            fprintf('Subject %s: Fold %d of %d.\n',subj,iFold,nFolds)
            [res.yPred{iFold},...
                res.yTest{iFold},...
                res.A{iFold},...
                res.intercept{iFold},...
                res.err{iFold}] = doRegFun(X,Y,folds(iFold,:),lambda,params.nCV,params.alpha,normType);
        end
        
        perf = mean(vertcat(res.err{:}));
        res.perf = perf;
        if doBinary
            yPred = vertcat(res.yPred{:});
            [~,~,~,res.AUC] = perfcurve(Y,yPred,true);
            %            res.AUC = compute_auc(yPred,Y);
            AUC = res.AUC;
        else
            mse = mean(vertcat(res.err{:}).^2);
            res.mse = mse;
            sse = sum(vertcat(res.err{:}).^2);
            res.sse = sse;
            mae = mean(abs(vertcat(res.err{:})));
            res.mae = mae;
            
            ssr = sum((vertcat(res.yTest{:})-vertcat(res.yPred{:})).^2);
            %                 sst = var(vertcat(res(t).yTest{:}))*length(vertcat(res(t).yTest{:}));
            sst = sum((vertcat(res.yTest{:}) - mean(vertcat(res.yTest{:}))).^2);
            
            r2  = 1 - ssr/sst;
            res.r2 = r2;
            
            r = corr(vertcat(res.yPred{:}),vertcat(res.yTest{:}));
            res.r = r;
            
        end
    end
    
    subject       = subj;
    params.lambda = lambda;
    if saveOutput
        save(fname,'res','Y','objLocs','params','perf','tal','AUC','mse','sse','mae','r2','r');
    end
catch e
    fname = fullfile(saveDir,[subj '_lasso_error.mat']);
    save(fname,'e')
end

function [stats,lambda] = calcLambda(X,Y,doBinary,nCV,alpha)
if isempty(nCV)
    nCV = round(length(Y)/2);
end

if doBinary
    
    [~,stats] = lassoglm(X,Y,'binomial','CV', nCV, 'NumLambda', 50,'alpha',alpha);
    
    % if the best model is a model with no features, then the data
    % are not all predictive of the behavior. I set a high lambda,
    % making it hard for any non-zero weights to get
    % through. Otherwise, you can end up with odd behavior on some
    % of the folds.
    if stats.DF(stats.Index1SE) == 0
        lambda = 5;
    else
        lambda = stats.Lambda1SE;
    end
    
else
    % center x
    xHat = mean(X, 1)';
    xCentered = X' - xHat * ones(1,size(X',2));
    
    % center y
    intercept = mean(Y);
    yCentered = round(Y - intercept,14);
    
    % compute optimal lamda
    [~, stats] = lasso(xCentered', yCentered, 'CV', nCV, 'NumLambda', 50,'alpha',alpha);
    lambda     = stats.LambdaMinMSE;
    
end

function [boxConstraint] = calcBoxConstraint(X,Y,nCV)
if isempty(nCV)
    nCV = round(length(Y)/2);
end

% z0 = [0,0];
% opts = optimset('TolX',1e-1,'TolFun',1e-1);
opts = optimset('TolX',5e-4,'TolFun',5e-4);
% fun = @(z)SVM_min_fn(X,Y,exp(z(1)),exp(z(2)),nCV);



m = 10;
fval = zeros(m,1);
z = zeros(m,1);
for j = 1:m;
    c = cvpartition(Y,'KFold',nCV);
    % minfn = @(z)kfoldLoss(fitcsvm(X,Y,'CVPartition',c,...
    % 'KernelFunction','linear','standardize',true,'BoxConstraint',exp(z(2)),...
    % 'KernelScale',exp(z(1))));
    minfn = @(z)kfoldLoss(fitcsvm(X,Y,'CVPartition',c,...
        'KernelFunction','linear','standardize',true,'BoxConstraint',exp(z)));
    [searchmin fval(j)] = fminsearch(minfn,randn(1,1),opts);
    z(j,:) = exp(searchmin)
end
z = z(fval == min(fval),:)
z = mean(z,1);

boxConstraint = z;

return
[searchmin fval] = fminsearch(minfn,randn(2,1),opts)

[z_opt,Crit] = fminsearch(fun,z0,opts);
[~, Features_opt] = fun(z_opt);

%************ Get optimal results **************
Acc = 1 - Crit;                       % Accuracy for model
rbf_sigma = exp(z_opt(1));
boxconstraint = exp(z_opt(2));
disp(sprintf('Max Acc: %2.2f, RBF sigma: %1.2f. Boxconstraint: %1.2f',Acc,rbf_sigma,boxconstraint))



function [Crit Features] = SVM_min_fn(X,Y,rbf_sigma,boxConstraint,nCV)
direction = 'backward';
opts = statset('display','iter');
kernel = 'rbf';

disp(sprintf('RBF sigma: %1.4f. Boxconstraint: %1.4f',rbf_sigma,boxConstraint))
c = cvpartition(Y,'k',nCV);
opts = statset('display','iter','TolFun',1e-3);
fun = @(x_train,y_train,x_test,y_test)SVM_class_fun(x_train,y_train,x_test,y_test,kernel,rbf_sigma,boxConstraint);
[fs,history] = sequentialfs(fun,X,Y,'cv',c,'direction',direction,'options',opts);

Features = find(fs==1);        % Features selected for given sigma and C
[Crit,h] = min(history.Crit);  % Mean classification error

function MCE = SVM_class_fun(x_train,y_train,x_test,y_test,kernel,rbf_sigma,boxconstraint)
svmStruct = svmtrain(x_train,y_train,'Kernel_Function','rbf','rbf_sigma',rbf_sigma,'boxconstraint',boxconstraint);
y_fit = svmclassify(svmStruct,x_test);
C = confusionmat(y_test,y_fit);
N = sum(sum(C));
MCE = N - sum(diag(C)); % No. misclassified sample


function penalty = calcPenalty(X,Y,nCV)
% Calculate the L2 penalty parameter. Following logic from Marc Coutanche's
% optimal_penalty_search, for find the best parameter across a broad range
% of values. Then focus in on the best value for a more precise measure.
if isempty(nCV)
    nCV = round(length(Y)/2);
end

% find best penatlty across a wide range
starting_vals = [0 0.01 0.1 1.0 10 100 1000 10000];
auc_pen = NaN(length(starting_vals),nCV);
for pVal = 1:length(starting_vals)
    thisPen = starting_vals(pVal);
    inds    = crossvalind('Kfold',size(X,1),nCV);
    for cv = 1:nCV
        
        % train data for this cv
        xTrain = [ones(sum(inds~=cv),1) X(inds~=cv,:)];
        yTrain = Y(inds~=cv);
        
        % test data for this cv
        xTest = [ones(sum(inds==cv),1) X(inds==cv,:)];
        yTest = Y(inds==cv);
        
        % train on this date with this cv
        out = logRegFun(yTrain',xTrain',thisPen);
        
        % test on held out
        p = exp(out.weights' * xTest')./(1+exp(out.weights' * xTest'));
        auc_pen(pVal,cv) = compute_auc(p,yTest);
    end
end

% pick best of the wide range
[maxBroad,bestInd] = max(nanmean(auc_pen,2));
bestBroadPen       = starting_vals(bestInd);

% search around the best value
if find(bestInd) == 0;
    pen_vals = linspace(0,starting_vals(find(bestInd)+1)/2,10);
elseif find(bestInd) == length(starting_vals)
    pen_vals = linspace(bestBroadPen/2,starting_vals(bestInd)*1.5,10);
else
    pen_vals = linspace(bestBroadPen/2,starting_vals(find(bestInd)+1)/2,10);
end

% loop over all new values
auc_pen_refined = NaN(length(pen_vals),nCV);
for pVal = 1:length(pen_vals)
    thisPen = pen_vals(pVal);
    inds    = crossvalind('Kfold',size(X,1),nCV);
    for cv = 1:nCV
        
        % train data for this cv
        xTrain = [ones(sum(inds~=cv),1) X(inds~=cv,:)];
        yTrain = Y(inds~=cv);
        
        % test data for this cv
        xTest = [ones(sum(inds==cv),1) X(inds==cv,:)];
        yTest = Y(inds==cv);
        
        % train on this date with this cv
        out = logRegFun(yTrain',xTrain',thisPen);
        
        % test on held out
        p = exp(out.weights' * xTest')./(1+exp(out.weights' * xTest'));
        auc_pen_refined(pVal,cv) = compute_auc(p,yTest);
    end
end

% pick best of the narrow range
[maxRefined,bestInd] = max(nanmean(auc_pen_refined,2));
bestRefinedPen       = pen_vals(bestInd);

% pick best of both
[~,ind] = max([maxBroad maxRefined]);
penalties = [bestBroadPen bestRefinedPen];
penalty = penalties(ind);


function [yPred,yTest,A,intercept,err] = doRegFun(X,Y,trainInds,lambda,nCV,alpha,normType)
%
% This does the regression.
% X = # trials x # features
% Y = # trials vector of responses
% trainInds = logical vector of training/test
%

if isempty(nCV)
    nCV = round(length(Y)/2);
end

doBinary = false;
if islogical(Y)
    doBinary = true;
end

% We will do binary classification Y is logical, which is currently the
% default
if doBinary
    
    % I'm under sampling the larger class so that we have equal numbers.
    % This isn't a great idea of more skewed dataset, but since we are
    % using a median threshold, it doesn't really matter here. Should also
    % do multiple rounds.
    yTrainBool = Y(trainInds);
    xTrain     = X(trainInds,:)';
    
    % if no lambda given, calculate lambda for this fold.
    if isempty(lambda)
        if strcmp(normType,'L1')
            [~, stats] = lassoglm(xTrain',yTrainBool,'binomial','CV', nCV, 'NumLambda', 25,'alpha',alpha);
            
            % if the best model is a model with no features, then the data
            % are not all predictive of the behavior. I set a high lambda,
            % making it hard for any non-zero weights to get
            % through. Otherwise, you can end up with odd behavior on some
            % of the folds.
            if stats.DF(stats.Index1SE) == 0
                lambda = 5;
            else
                lambda = stats.Lambda1SE;
            end
        elseif strcmp(normType,'L2')
            lambda = calcPenalty(xTrain',yTrainBool,nCV);
        elseif strcmpi(normType,'SVM')
            lambda = calcBoxConstraint(xTrain',yTrainBool,nCV);
        end
    end
    
    % figure out which observations to remove from larger class
    numToRemove = sum(yTrainBool) - sum(~yTrainBool);
    
    nSubSamps = 1;
    if numToRemove ~= 0
        nSubSamps = 11;
    end
    
    yPreds     = NaN(sum(~trainInds),nSubSamps);
    As         = NaN(size(X,2),nSubSamps);
    intercepts = NaN(1,nSubSamps);
    
    for nSub = 1:nSubSamps
        
        toRemove = [];
        yTrainBool = Y(trainInds);
        if numToRemove > 0
            toRemove = randsample(find(yTrainBool),abs(numToRemove));
        elseif numToRemove < 0
            toRemove = randsample(find(~yTrainBool),abs(numToRemove));
        end
        
        % training set x
        xTrain = X(trainInds,:)';
        xTrain(:,toRemove) = [];
        
        % training set y
        yTrainBool(toRemove) = [];
        
        link = 'logit';
        if strcmp(normType,'L1')
            [As(:,nSub), stats] = lassoglm(xTrain',yTrainBool,'binomial','Lambda',lambda,'alpha',alpha);
            intercepts(nSub)    = stats.Intercept;
        elseif strcmpi(normType','L2')
            xTrain           = [ones(sum(size(xTrain,2)),1) xTrain'];
            out              = logRegFun(yTrainBool',xTrain',lambda);
            intercepts(nSub) = out.weights(1);
            As(:,nSub)       = out.weights(2:end);
        elseif strcmpi(normType,'svm')
            link = 'identity';
            % SVMModel = fitcsvm(xTrain',yTrainBool,'KernelFunction','linear','KernelScale',kScale,'BoxConstraint',lambda,'standardize',true);
            SVMModel = fitcsvm(xTrain',yTrainBool,'BoxConstraint',lambda);
            %'standardize',true);
            % [ScoreSVMModel,Info] = fitPosterior(SVMModel);
            intercepts(nSub) = SVMModel.Bias;
            As(:,nSub)       = SVMModel.Beta;
            % slopes = Info.Slope;
            % B1 = [ScoreSVMModel.Bias+Info.Intercept;ScoreSVMModel.Beta];
            % keyboard
        end
        
        %
        % scores=(xTest')*ScoreSVMModel.Beta + ScoreSVMModel.Bias
        %         1./(1+ exp(Info.Slope*scores + Info.Intercept))
        %
        %
        % 1./(1+ exp(Info.Slope*scores + Info.Intercept))
        
        % testing set
        xTest = X(~trainInds,:)';
        yTest = Y(~trainInds);
        
        % predict
        B1 = [intercepts(nSub);As(:,nSub)];
        yPred = glmval(B1,xTest',link);
        
        % In cases where the model is empty, I'm going to randomly move the
        % prediction off of the mean in either the positive or negative
        % direction. This is a little weird, but I don't want it to always say
        % that when the model is empty, choose recalled.
        %if nnz(As(:,nSub)) == 0
        %    yPred = yPred + (randsample([-1 1],length(yPred),true) * .0001)';
        %end
        
        yPreds(:,nSub) = yPred;
    end
    
    yPred     = mean(yPreds,2);
    A         = mean(As,2);
    intercept = mean(intercepts);
    
    % see if the predictions match the actual results
    if ~strcmpi(normType,'svm')
        err = (yPred > mean(yTrainBool)) == Y(~trainInds);
    else
        err = (yPred > 0) == Y(~trainInds);
    end
    
    %     err = mean((yPreds > mean(yTrainBool)) == repmat(Y(~trainInds),1,nSub),2);
    
    % if Y is not logical, do regression
else
    
    % training set x
    xTrain = X(trainInds,:)';
    %     xHat = mean(xTrain, 2);
    %     xTrain = xTrain - xHat * ones(1,size(xTrain,2));
    
    % training set y
    yTrain = Y(trainInds);
    intercept = mean(yTrain);
    intercept = mean(Y);
    %     yTrain = round(yTrain - intercept,14);
    
    % compute model
    if isempty(lambda)
        [A_lasso, stats] = lasso(xTrain', yTrain, 'CV', nCV, 'NumLambda', 50, 'alpha',alpha);
        A = A_lasso(:,stats.IndexMinMSE);
    else
        A = lasso(xTrain', yTrain, 'Lambda',lambda, 'alpha',alpha);
    end
    
    % testing set
    xTest = X(~trainInds,:)';
    yTest = Y(~trainInds);
    
    % Double check this
    %     yPred = (xTest - xHat*ones(1,sum(~trainInds)))' * A + intercept;
    yPred = (xTest)' * A + intercept;
    %     err = mean((yTest - yPred).^2);
    err = yTest - yPred;
    
    
    
    
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
        fname = fullfile(subjPath,'RAM_YC1_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
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
        tInds = powParams.timeBins(:,1) >= timeBins(t,1) & powParams.timeBins(:,2) <= timeBins(t,2);
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














