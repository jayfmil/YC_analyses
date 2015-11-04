function [perf,AUC,subject,params] = YC1_runMulti_subj_liblinear(subj,params,saveDir)
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
    
    %%% THIS MIGHT BE WRONG %%%
    % and filter to encoding period of interest
    [~,firstIdx,~] = unique([session(eventsToUse)' [events(eventsToUse).blocknum]'],'rows','first');
    if any(strcmp({'first','second'},params.encPeriod))        
        if strcmp(params.encPeriod,'first')
            eventsToUse(firstIdx+1) = false;
        elseif strcmp(params.encPeriod,'second')
            eventsToUse(firstIdx) = false;
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
        [powerData,sessInds] = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params);
        powerData = permute(powerData,[3 1 2 4]);
        if params.savePower
            powFile = fullfile(powDir,[subj '_binnedPower.mat']);
            save(powFile,'powerData','params','sessInds')
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
                    fprintf('Subject %s: Computing optimal c using %s regularization.\n',subj,normType)
                    lambda = calcPenalty(X,Y,params.nCV,normType,sessInds,folds);
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
                    res(t).err{iFold}] = doRegFun(X,Y,folds,iFold,lambda,normType,sessInds);
                
            end
            
            
            perf(t) = mean(vertcat(res(t).err{:}));
            res(t).perf = perf(t);
            yPred = vertcat(res(t).yPred{:});
            [~,~,~,res(t).AUC] = perfcurve(Y,yPred,true);
            %                res(t).AUC = compute_auc(yPred,Y);
            AUC(t) = res(t).AUC;
            
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
                res.err{iFold}] = doRegFun(X,Y,folds(iFold,:),lambda,params.alpha,normType);
        end
        
        perf = mean(vertcat(res.err{:}));
        res.perf = perf;
        
            yPred = vertcat(res.yPred{:});
            [~,~,~,res.AUC] = perfcurve(Y,yPred,true);
            %            res.AUC = compute_auc(yPred,Y);
            AUC = res.AUC;
      
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

function penalty = calcPenalty(X,Y,folds,sessInds,normType)
% Calculate the penalty parameter c.

Y = double(Y);
Y(Y==0) = -1;
C       = 2.^(-20:2:20);
nCV     = size(folds,1);
auc_pen = NaN(length(C),1);

for pVal = 1:length(C)
    thisPen = C(pVal);    

    % hold test output
    dec_values = [];
    labels = [];
    
    % set parameters for either L1 or L2 with the current c
    if strcmpi(normType,'L1')
        liblin_param = ['-c ' num2str(thisPen) ' -s 6'];
    elseif strcmpi(normType,'L2')
        liblin_param = ['-c ' num2str(thisPen) ' -s 0'];
    end
    
    for cv = 1:nCV
        inds = folds(cv,:);
        
        % train data for this cv
        xTrain = [ones(sum(inds),1) X(inds,:)];
        yTrain = Y(inds);        
        [xTrain,m,s] = standardize(xTrain,sessInds(inds));
        
        % test data for this cv
        xTest = [ones(sum(~inds),1) X(~inds,:)];
        yTest = Y(~inds);
        xTest = standardize_test(xTest,sessInds(~inds),m,s);
        
        % train on this data with this cv        
        model = train(double(yTrain),sparse(xTrain),liblin_param);
        
        % predict with model
        [pred, acc, dec] = predict(double(yTest),sparse(xTest),model);
        
        % keep of track of output over all folds
        if model.Label(1) < 0;
            dec = dec * -1;
        end
        dec_values = vertcat(dec_values, dec);
        labels = vertcat(labels, yTest);
        
    end    
    
    % find AUC for this value of C
    [~,~,~,auc_pen(pVal)] = perfcurve(labels,dec_values,1);    
end

% return C with highest AUC
penalty = C(find(auc_pen==max(auc_pen),1,'first'));

function [yProbs,yTest,A,err] = doRegFun(X,Y,folds,iFold,lambda,normType,sessInds)
%X,Y,folds,iFold,lambda,params.nCV,normType
% This does the regression.
% X = # trials x # features
% Y = # trials vector of responses
% trainInds = logical vector of training/test
%

% get data for this fold
trainInds  = folds(iFold,:);
sessions   = sessInds(trainInds);
subFolds   = folds(setdiff(1:size(folds,1),iFold),trainInds);
yTrainBool = Y(trainInds);
xTrain     = X(trainInds,:);

% if no lambda given, calculate lambda for this fold.
if isempty(lambda)    
    lambda = calcPenalty(xTrain,yTrainBool,subFolds,sessions,normType);         
end

% figure out which observations to remove from larger class to make the
% sizes equal
numToRemove = sum(yTrainBool) - sum(~yTrainBool);
nSubSamps = 1;
if numToRemove ~= 0
    nSubSamps = 11;
end

% subsampling
yProbs     = NaN(sum(~trainInds),nSubSamps);
yPreds       = NaN(sum(~trainInds),nSubSamps);
As         = NaN(size(X,2)+1,nSubSamps);
intercepts = NaN(1,nSubSamps);

% loop over each sub sample
for nSub = 1:nSubSamps
    
    % remove random training
    toRemove = [];
    yTrainBool = Y(trainInds);
    if numToRemove > 0
        toRemove = randsample(find(yTrainBool),abs(numToRemove));
    elseif numToRemove < 0
        toRemove = randsample(find(~yTrainBool),abs(numToRemove));
    end
    
    % training set x
    sessIndsTmp = sessions;
    sessIndsTmp(toRemove) = [];
    xTrain = X(trainInds,:);
    xTrain(toRemove,:) = [];
    
    % training set y
    yTrainBool(toRemove) = [];
    yTrain = double(yTrainBool);
    yTrain(yTrain==0) = -1;
    
    % standardize training data
    xTrain = [ones(size(xTrain,1),1) xTrain];      
    [xTrain,m,s] = standardize(xTrain,sessIndsTmp);
        
    % set parameters for either L1 or L2
    if strcmpi(normType,'L1')
        liblin_param = ['-c ' num2str(lambda) ' -s 6'];
    elseif strcmpi(normType,'L2')
        liblin_param = ['-c ' num2str(lambda) ' -s 0'];
    end    
    
    % train model
    model = train(double(yTrain),sparse(xTrain),liblin_param);
    As(:,nSub) = model.w;    
    
    % testing set
    xTest = X(~trainInds,:);
    yTest = double(Y(~trainInds));
    yTest(yTest==0) = -1;
    
    % standardize testing data
    xTest = [ones(size(xTest,1),1) xTest];    
    xTest = standardize_test(xTest,sessInds(~trainInds),m,s);
    
    % predict
    [pred, acc, yProb] = predict(double(yTest),sparse(xTest),model,'-b 1');
        
    yProbs(:,nSub) = yProb(:,1);
    yPreds(:,nSub) = pred;
    
end

yProbs     = mean(yProbs,2); 
A          = mean(As,2);
% intercept = mean(intercepts);

err = mean(yPreds == repmat(yTest,1,nSub),2);


function [X,m,s] = standardize(X,sessInds)

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
X(:,1) = 1;


function [X] = standardize_test(X,sessInds,m,s)

sessions = unique(sessInds);
for sess = 1:length(sessions)
    
    X(sessInds==sessions(sess),:) = bsxfun(@minus, X(sessInds==sessions(sess),:),...
        m(sess,:));
    
    X(sessInds==sessions(sess),:) = bsxfun(@rdivide, X(sessInds==sessions(sess),:),...
        s(sess,:));        
end
X(:,1) = 1;



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
        subjPow = cat(3,subjPow,pow);
        sessInds = cat(1,sessInds,ones(nE,1)*s);
    end
    
    if length(eventsToUse) ~= size(subjPow,3)
        fprintf('Number of events does not match size of power matrix for %s!.\n',subj)
        return
    end
    subjPow  = subjPow(:,:,eventsToUse);
    sessInds = sessInds(eventsToUse);    
        
    
%     % average frequencies
%     if nFreqs ~= length(powParams.params.pow.freqs)
%         tmpPower = NaN(nFreqs,size(subjPow,2),size(subjPow,3));
%         for f = 1:nFreqs
%             fInds = powParams.params.pow.freqs >= freqBins(f,1) & powParams.params.pow.freqs < freqBins(f,2);
%             tmpPower(f,:,:) = nanmean(subjPow(fInds,:,:),1);
%         end
%         subjPow = tmpPower;
%     end
    
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














