function [perf,AUC,subject,params] = YC1_runMulti_subj(subj,params,saveDir)
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
        powerData = powerData.powerData;
    else
        powerData = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params);
        powerData = permute(powerData,[3 1 2 4]);
        if params.savePower
            powFile = fullfile(powDir,[subj '_binnedPower.mat']);
            save(powFile,'powerData','params')
        end
    end
    
    % size of feature matrix
    nElecs = size(powerData,4);
    nTimes = size(powerData,3);
    nFreqs = size(powerData,2);
    
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
    % Note: this is not influenced by params.nCV. params.nCV is the number
    % of cross validation folds to estimate lambda with lassoglm. For the
    % training test iterations, we are always doing leave one object out.
    [trials,~,trialInds] = unique([session(eventsToUse)' [events(eventsToUse).blocknum]'],'rows');
    nFolds = size(trials,1);
    folds = false(nFolds,size(trialInds,1));
    for iFold = 1:nFolds
        folds(iFold,:) = trialInds ~= iFold;
    end
    
    % response data
    Y = [events(eventsToUse).testError]';
    if doBinary
        Y  = Y < median(Y);
    end
        
    % use hagai's error metric instead of the normal one?
    if isfield(params,'useHagai') && params.useHagai == 1
        clusteredLabels = load('/home1/jfm2/matlab/YC_analyses/Multivariate/clusteredLabels2.mat');
        sInd = strcmp(clusteredLabels.subjects,subj);
        hagaiLabel = clusteredLabels.binaryBehaviorMetric{sInd};
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
        Ytmp = reshape(Y,2,[]);
        Y = reshape(Ytmp(:,randOrder),[],1);
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
            X = reshape(squeeze(powerData(:,:,t,:)),size(powerData,1),nFreqs*nElecs);
            
            % see if we are precomputing lambda based on all the data
            lambda = [];
            if params.crossValStrictness == 0
                if ~isempty(params.lambda)
                    lambda = params.lambda(t);
                else
                    fprintf('Subject %s: Computing optimal lambda.\n',subj)
                    [stats,lambda] = calcLambda(X,Y,doBinary,params.nCV);
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
                    res(t).err{iFold}] = lassoReg(X,Y,folds(iFold,:),lambda,params.nCV);
            end  
            perf(t) = mean(vertcat(res(t).err{:}));
            res(t).perf = perf(t);
            if doBinary     
                yPred = vertcat(res(t).yPred{:});                
%                 [res(t).xAUC,res(t).yAUC,~,res(t).AUC,res(t).optPoint] = perfcurve(Y,yPred,true);
                res(t).AUC = compute_auc(yPred,Y);
                AUC(t) = res(t).AUC;                
            end
        end
        lambda = lambdas;
        
    % if using all time points in one model, current what it is set to do
    else
        
        % reshape into # trials x # features
        X = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs);
        
        % see if we are precomputing lambda based on all the data
        lambda = [];
        if params.crossValStrictness == 0
            if ~isempty(params.lambda)
                lambda = params.lambda;
            else
                fprintf('Subject %s: Computing optimal lambda.\n',subj)
                [stats,lambda] = calcLambda(X,Y,doBinary,params.nCV);
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
                res.err{iFold}] = lassoReg(X,Y,folds(iFold,:),lambda,params.nCV);
        end
        perf = mean(vertcat(res.err{:}));
        res.perf = perf;
        if doBinary
            yPred = vertcat(res.yPred{:});
%             [res.xAUC,res.yAUC,~,res.AUC,res.optPoint] = perfcurve(Y,yPred,true);
            res.AUC = compute_auc(yPred,Y);
            AUC = res.AUC;            
        end
    end
        
    subject       = subj;
    params.lambda = lambda;
    if saveOutput        
        save(fname,'res','Y','objLocs','params','perf','tal','AUC');
    end
catch e
    fname = fullfile(saveDir,[subj '_lasso_error.mat']);
    save(fname,'e')
end

function [stats,lambda] = calcLambda(X,Y,doBinary,nCV)
if isempty(nCV)
    nCV = round(length(Y)/2);
end

if doBinary
  
    [~,stats] = lassoglm(X,Y,'binomial','CV', nCV, 'NumLambda', 50);
    lambda    = stats.LambdaMinDeviance;
    
else
    % center x
    xHat = mean(X, 1)';
    xCentered = X' - xHat * ones(1,size(X',2));
    
    % center y
    intercept = mean(Y);
    yCentered = round(Y - intercept,14);
    
    % compute optimal lamda
    [~, stats] = lasso(xCentered', yCentered, 'CV', nCV, 'NumLambda', 50);
    lambda     = stats.LambdaMinMSE;
end

function [yPred,yTest,A,intercept,err] = lassoReg(X,Y,trainInds,lambda,nCV)
%
% This does lasso. 
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
    % using a median threshold, it doesn't really matter here
    yTrainBool = Y(trainInds);
    
    % figure out which observations to remove from larger class
    numToRemove = sum(yTrainBool) - sum(~yTrainBool);
    toRemove = [];
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
    
    % compute model
    if isempty(lambda)
        [A_lasso, stats] = lassoglm(xTrain',yTrainBool,'binomial','CV', nCV, 'NumLambda', 50);
        
        % get the best cooefficients and intercept
        A = A_lasso(:,stats.IndexMinDeviance);
        intercept = stats.Intercept(stats.IndexMinDeviance);
    else
        [A, stats] = lassoglm(xTrain',yTrainBool,'binomial','Lambda',lambda);
        intercept = stats.Intercept;
    end
    
    % testing set
    xTest = X(~trainInds,:)';
    yTest = Y(~trainInds);
    
    % predict
    B1 = [intercept;A];
    yPred = glmval(B1,xTest','logit');
    
    % see if the predictions match the actual results
    err = (yPred > mean(yTrainBool)) == Y(~trainInds);
    
% if Y is not logical, do regression
else
    
    % training set x
    xTrain = X(trainInds,:)';    
    xHat = mean(xTrain, 2);
    xTrain = xTrain - xHat * ones(1,size(xTrain,2));
    
    % training set y
    yTrain = Y(trainInds);
    intercept = mean(yTrain);
    yTrain = round(yTrain - intercept,14);
    
    % compute model
    if isempty(lambda)
        [A_lasso, stats] = lasso(xTrain', yTrain, 'CV', nCV, 'NumLambda', 50);
        A = A_lasso(:,stats.IndexMinMSE);
    else
        A = lasso(xTrain', yTrain, 'Lambda',lambda);
    end
    
    % testing set
    xTest = X(~trainInds,:)';
    yTest = Y(~trainInds);
    
    % Double check this
    yPred = (xTest - xHat*ones(1,sum(~trainInds)))' * A + intercept;
    err = mean((yTest - yPred).^2);
           
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














