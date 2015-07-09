function [perf,subject,params] = YC1_runMulti_subj(subj,params,saveDir)

perf = [];
subject = [];

try
    
    % load subject electrode locations
    tal = getBipolarSubjElecs(subj,1,1,1);
    if ~isfield(tal,'locTag') || ~any(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
        fprintf('No MTL electrodes for %s.\n',subj)
        return
    end    
    
    % load events so we can filter into our conditions of interest
    config = RAM_config('RAM_YC1');
    
    % Setting time bins for convenience:
    tEnds = (config.distributedParams.timeWin:...
        config.distributedParams.timeStep:...
        config.postMS+config.priorMS)-config.priorMS;
    tStarts = tEnds - config.distributedParams.timeWin + 1;
    config.distributedParams.timeBins = [tStarts' tEnds'];
    
    % load events
    [events] = RAM_loadEvents(subj,[],'RAM_YC1','events', config);
    
    % add the test error to the learning trials
    events = addErrorField(events);
    
    % convert to the session to an double from a string
    session = NaN(1,length(events));
    for e = 1:length(session)
        session(e) = str2double(events(e).session);
    end
    
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
    
    % see if we are passing in a lambda value. If so, will not compute
    % optimal one unless params.crossValStrictness is true
    if isfi
    
    % load power for all electrodes
    if ~params.useCorrectedPower
        powerData = loadAllPower(tal,subj,events,freqBins,timeBins,config,eventsToUse);
    else
        powerData = loadResids_locs(tal,subj,freqBins,timeBins,config,eventsToUse);
    end
    powerData = permute(powerData,[3 1 2 4]);
    if params.savePower
        powDir = fullfile(saveDir,'power');
        if ~exist(powDir,'dir')
            mkdir(powDir);
        end
        powFile = fullfile(powDir,[subj '_binnedPower.mat']);
        save(powFile,'powerData','params')
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
    
    if doPermute
        randOrder = randperm(size(trials,1));
        randOrder = repmat(randOrder,2,1);
        randOrder = randOrder(:);
        Y = Y(randOrder);
    end
    
    objLocs = vertcat(events(eventsToUse).objLocs);
        
    % We can model time points seperately, so # features = # freqs x # elecs,
    % or we can model it all together, so # features = # times x # freqs x #
    % elecs.    
    res = [];
    if modelEachTime
        perf = NaN(1,nTimes);
        for t = 1:nTimes
            
            % reshape into # trials x # features    
            X = reshape(squeeze(powerData(:,:,t,:)),size(powerData,1),nFreqs*nElecs);
            
            % see if we are precomputing lambda based on all the data
            lambda = [];
            if params.crossValStrictness == 0
                if ~isempty(params.lambda)
                    lambda = params.lambda;
                else
                    [stats,lambda] = calcLambda(X,Y,doBinary);
                end
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
                    res(t).err{iFold}] = lassoReg(X,Y,folds(iFold,:),lambda);
            end  
            perf(t) = mean(vertcat(res(t).err{:}));
        end
        
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
                [stats,lambda] = calcLambda(X,Y,doBinary);
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
                res.err{iFold}] = lassoReg(X,Y,folds(iFold,:),lambda);
        end
        perf = mean(vertcat(res.err{:}));
    end
      
    subject       = subj;
    params.lambda = lambda;
    if saveOutput
        fname = fullfile(saveDir,[subj '_lasso.mat']);
        save(fname,'res','Y','objLocs','params','perf');
    end
catch e
    fname = fullfile(saveDir,[subj '_error_lasso.mat']);
    save(fname,'e')
end

function [stats,lambda] = calcLambda(X,Y,doBinary)
if doBinary
  
    [~,stats] = lassoglm(X,Y,'binomial','CV', round(length(Y)/2), 'NumLambda', 50);
    lambda    = stats.LambdaMinDeviance;
    
else
    % center x
    xHat = mean(X, 1)';
    xCentered = X' - xHat * ones(1,size(X',2));
    
    % center y
    intercept = mean(Y);
    yCentered = round(Y - intercept,14);
    
    % compute optimal lamda
    [~, stats] = lasso(xCentered', yCentered, 'CV', round(length(Y)/2), 'NumLambda', 50);
    lambda     = stats.LambdaMinMSE;
end

function [yPred,yTest,A,intercept,err] = lassoReg(X,Y,trainInds,lambda)
%
% This does lasso. 
% X = # trials x # features
% Y = # trials vector of responses
% trainInds = logical vector of training/test
%

doBinary = false;
if islogical(Y)
    doBinary = true;
end

% We will do binary classification Y is logical, which is currently the
% default
if doBinary
    
    % I'm under sampling the larger class so that we have equal numbers.
    % I'm doing this so that I can know that .5 will be my cutoff for
    % classification. Maybe there is a something else I should be doing?
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
        [A_lasso, stats] = lassoglm(xTrain',yTrainBool,'binomial','CV', 5, 'NumLambda', 10);
        
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
    err = round(yPred) == Y(~trainInds);
    
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
        [A_lasso, stats] = lasso(xTrain', yTrain, 'CV', 5, 'NumLambda', 25);
        A = A_lasso(:,stats.IndexMinMSE);
    else
        A = lasso(xTrain', yTrain, 'Lambda',lambda);
    end
    
    % testing set
    xTest = X(~trainInds,:)';
    yTest = Y(~trainInds);
    
    % I'm predicting using the intercept we computed above, not the one
    % returned from lasso. Is that wrong?
    yPred = (xTest - xHat*ones(1,sum(~trainInds)))' * A + intercept;
    err = mean((yTest - yPred).^2);
           
end



function powerData = loadAllPower(tal,subj,events,freqBins,timeBins,config,eventsToUse)

nFreqs = size(freqBins,1);
nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
nElecs = length(tal);
powerData = NaN(nFreqs,nTimes,nEvents,nElecs);

for e = 1:nElecs
    elecNum = tal(e).channel;
    
    % load power for all sessions. Power should aleady have been
    % created or else error
    [distOut] = RAM_dist_func(subj,[],elecNum,'RAM_YC1','events', 0, ...
        @doNothing, config.distributedFunctionLabel, config.distributedParams, 1, 1, events);
    
    sessInds = distOut.sessInds;
    subjMean = distOut.meanBasePow;
    subjStd = distOut.stdBasePow;
    subjPow = distOut.pow;
    
    % a couple sessions weren't zscored, so do it here. I should double
    % check that this is right
    if isempty(subjStd)
        fprintf('power not zscored for %s\n',subj)
        
        sI = unique(sessInds);
        zpow = NaN(size(subjPow));
        for s = 1:length(sI)
            inds = sessInds == sI(s);
            subjMean = nanmean(squeeze(nanmean(subjPow(:,:,inds),2)),2);
            subjMean = repmat(subjMean,[1 size(subjPow,2), size(subjPow,3)]);
            subjStd = nanstd(squeeze(nanmean(subjPow(:,:,inds),2)),[],2);
            subjStd = repmat(subjStd,[1 size(subjPow,2), size(subjPow,3)]);
            zpow(:,:,inds) = (subjPow(:,:,inds) - subjMean).*subjStd;
        end
        subjPow = zpow;
    end
    subjPow = subjPow(:,:,eventsToUse);
    
    % average frequencies
    tmpPower = NaN(nFreqs,size(subjPow,2),size(subjPow,3));
    for f = 1:nFreqs
        fInds = config.distributedParams.freQ >= freqBins(f,1) & config.distributedParams.freQ < freqBins(f,2);
        tmpPower(f,:,:) = nanmean(subjPow(fInds,:,:),1);
    end
    subjPow = tmpPower;
    
    % average times
    tmpPower = NaN(nFreqs,nTimes,size(subjPow,3));
    for t = 1:nTimes
        tInds = config.distributedParams.timeBins(:,1) >= timeBins(t,1) & config.distributedParams.timeBins(:,2) < timeBins(t,2);
        tmpPower(:,t,:) = nanmean(subjPow(:,tInds,:),2);
    end
    powerData(:,:,:,e) = tmpPower;
end

function [powerData] = loadResids_locs(tal,subj,freqBins,timeBins,config,eventsToUse)

nFreqs = size(freqBins,1);
nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
nElecs = length(tal);
powerData = NaN(nFreqs,nTimes,nEvents,nElecs);

basePath  = '/data10/scratch/jfm2/YC1/multi/power/regress/';
subjPath  = fullfile(basePath,subj);

for e = 1:nElecs
    elecNum   = tal(e).channel;
    fname     = sprintf('%s_elec_%d-%d_residuals.mat',subj,elecNum(1),elecNum(2));
    
    if ~exist(fullfile(subjPath,fname),'file')
        error('Residuals file %s not found.\n',fname)
    else
        elecData = load(fullfile(subjPath,fname));
        resids   = elecData.resid;
        resids   = permute(resids,[3 2 1]);
        if size(resids,3) ~= sum(eventsToUse)
            keyboard
        end
        
        % average frequencies
        tmpPower = NaN(nFreqs,size(resids,2),size(resids,3));
        for f = 1:nFreqs
            fInds = config.distributedParams.freQ >= freqBins(f,1) & config.distributedParams.freQ < freqBins(f,2);
            tmpPower(f,:,:) = nanmean(resids(fInds,:,:),1);
        end
        resids = tmpPower;
        
        % average times
        tmpPower = NaN(nFreqs,nTimes,size(resids,3));
        for t = 1:nTimes
            tInds = config.distributedParams.timeBins(:,1) >= timeBins(t,1) & config.distributedParams.timeBins(:,2) < timeBins(t,2);
            tmpPower(:,t,:) = nanmean(resids(:,tInds,:),2);
        end
        powerData(:,:,:,e) = tmpPower;
        
    end
end



% reshape into number of observations (events) x number of features
% powerData = permute(powerData,[3 1 2 4]);
% powerData = reshape(powerData,size(powerData,1),nFreqs*nTimes*nElecs);


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














