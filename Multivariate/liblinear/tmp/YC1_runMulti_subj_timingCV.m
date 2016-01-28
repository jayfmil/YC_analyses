function [perf,AUC,subject,params] = YC1_runMulti_subj_timingCV(subj,params,saveDir)
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
    cvField       = params.cvField;   
    nestedCvField = params.nestedCvField;   
    prctileThresh = params.auc_prctileThresh;
    
    % load events
    events = get_sub_events('RAM_YC1',subj);
    
    % add the test error to the learning trials
    events  = addExtraFields(events);
    session = [events.session];
    
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
    
    %%% THIS MIGHT BE WRONG I SHOULD FIX IT %%%
    % and filter to encoding period of interest
    [~,firstIdx,~] = unique([session(eventsToUse)' [events(eventsToUse).blocknum]'],'rows','first');
    if any(strcmp({'first','second'},params.encPeriod))        
        if strcmp(params.encPeriod,'first')
            eventsToUse(firstIdx+1) = false;
        elseif strcmp(params.encPeriod,'second')
            eventsToUse(firstIdx) = false;
        end
    end    
    
    
    % load power for all electrodes
    powDir = fullfile(saveDir,'power');
    powFile = fullfile(powDir,[subj '_binnedPower.mat']);
    if ~exist(powDir,'dir')
        mkdir(powDir)
    end
    if params.loadPower
        powerData = load(powFile);
        timeLabel = powerData.timeLabel;
        sessInds  = powerData.sessInds;
        powerData = powerData.powerData;        
    else
        [powerData,sessInds] = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params);
        powerData = permute(powerData,[3 1 2 4]);
        timeLabel = NaN(1,size(powerData,2),size(powerData,3),size(powerData,4));
        for t = 1:size(timeLabel,3)
            timeLabel(:,:,t,:) = t;
        end
        if params.savePower
            powFile = fullfile(powDir,[subj '_binnedPower.mat']);
            save(powFile,'powerData','params','sessInds','timeLabel')
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
    [nFolds,folds,trials] = createFolds(session(eventsToUse),[events(eventsToUse).(cvField)]);
%     if strcmpi(params.encPeriod,'combined') || strcmpi(params.encPeriod,'average')
%         trialInds = trialInds(firstIdx);    
%         Y = Y(firstIdx);
%     end    

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
          
    % Unlike previous versions, if modelEachTime is true, this will not
    % create nTimes number of models. Instead, the model for fold each will
    % be chosen based on inner fold cross validation of both penalty and
    % time point.
    res = [];
    X   = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs*nItemObs);
    T   = reshape(squeeze(timeLabel),size(timeLabel,1),nFreqs*nTimes*nElecs*nItemObs);
    
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
        res.lambda,res.tBest,res.Cs,res.Ts,res.aucs] = deal(cell(nFolds,1));
    for iFold = 1:nFolds
        fprintf('Subject %s: Fold %d of %d.\n',subj,iFold,nFolds)
        [res.yProb{iFold},...
            res.yPred{iFold},...
            res.yTest{iFold},...
            res.A{iFold},...
            res.err{iFold},...
            res.lambda{iFold},...
            res.tBest{iFold},...
            res.Cs{iFold},...
            res.Ts{iFold},...
            res.aucs{iFold}] = doRegFun(X,Y,T,folds,iFold,lambda,normType,sessInds,events(eventsToUse),nestedCvField,prctileThresh,modelEachTime);
        fprintf('Subject %s: Fold %d of %d, %.3f running percent correct.\n',subj,iFold,nFolds,mean(vertcat(res.err{:})))
    end
    
    yProb = vertcat(res.yProb{:});
    [~,~,~,res.AUC] = perfcurve(Y,yProb,true);
    AUC = res.AUC;
    perf = mean(vertcat(res.err{:}));
    
    % Finally, create a model based on the full dataset for the patient. If
    % modelEachTime, do it for the optimal time
%     bestTime = find(AUC == max(AUC),1,'first');
%     X = reshape(squeeze(powerData(:,:,bestTime,:,:)),size(powerData,1),nFreqs*nElecs*nItemObs);
    resFull = doRegFullData(X,Y,T,folds,normType,sessInds,prctileThresh);
    
    
%     if modelEachTime
%         perf = NaN(1,nTimes);
%         lambdas = NaN(1,nTimes);
%         
% 
%         
%         for t = 1:nTimes
%             
%             % reshape into # trials x # features
%             X = reshape(squeeze(powerData(:,:,t,:,:)),size(powerData,1),nFreqs*nElecs*nItemObs);
%             
%             % see if we are precomputing lambda based on all the data. This
%             % happens if crossValStrictness is not 1. This is technically
%             % not the correct to do things
%             lambda = [];
%             if params.crossValStrictness == 0
%                 if ~isempty(params.lambda)
%                     lambda = params.lambda(t);
%                 else                    
%                     fprintf('Subject %s: Computing optimal c using %s regularization.\n',subj,normType)
%                     lambda = calcPenalty(X,Y,T,folds,sessInds,normType);
%                 end
%                 lambdas(t) = lambda;
%             end
%             
%             % will hold results from each fold
%             % yPred  = predicted class
%             % yTest  = actual class
%             % err    = correct or not
%             % lambda = value of penatly paramater for the fold
%             % Cs     = all possible penalty parameters tried
%             % aucs   = aucs for each penalty param when choosing the best C
%             [res(t).yPred,res(t).yTest,...
%              res(t).A,res(t).err,...
%              res(t).lambda,res(t).Cs,res(t).aucs] = deal(cell(nFolds,1));
%             
%             % run for each fold
%             for iFold = 1:nFolds
%                 fprintf('Subject %s: Time %d of %d, Fold %d of %d.\n',subj,t,nTimes,iFold,nFolds)
%                 [res(t).yPred{iFold},...
%                     res(t).yTest{iFold},...
%                     res(t).A{iFold},...
%                     res(t).err{iFold},...
%                     res(t).lambda{iFold},...
%                     res(t).Cs{iFold},...
%                     res(t).aucs{iFold}] = doRegFun(X,Y,folds,iFold,lambda,normType,sessInds);                
%             end             
%             
%             % calculate performance metrics
%             % percent correct
%             perf(t) = mean(vertcat(res(t).err{:}));
%             res(t).perf = perf(t);
%             
%             % auc
%             yPred = vertcat(res(t).yPred{:});
%             [~,~,~,res(t).AUC] = perfcurve(Y,yPred,true);
%             AUC(t) = res(t).AUC;            
%             
%             % save out the mean lambda for convenience?
%             lambdas(t) = mean([res(t).lambda{:}]);
%             
%         end
%         lambda = lambdas;
%         
%         % Finally, create a model based on the full dataset for the patient. If
%         % modelEachTime, do it for the optimal time       
%         bestTime = find(AUC == max(AUC),1,'first');
%         X = reshape(squeeze(powerData(:,:,bestTime,:,:)),size(powerData,1),nFreqs*nElecs*nItemObs);        
%         resFull = doRegFullData(X,Y,folds,normType,sessInds);
%         
%     % if using all time points in one model, current what it is set to do
%     else
%         
%         % reshape into # trials x # features
%         X = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs*nItemObs);
%         
%         % see if we are precomputing lambda based on all the data. Again,
%         % we shouldn't really do this
%         lambda = [];
%         if params.crossValStrictness == 0
%             if ~isempty(params.lambda)
%                 lambda = params.lambda;
%             else
%                 fprintf('Subject %s: Computing optimal c using %s regularization.\n',subj,normType)
%                 lambda = calcPenalty(X,Y,folds,sessInds,normType);
%             end
%         end        
%         
%         % run for each fold
%         [res.yPred,res.yTest,res.A,res.intercept,res.err,...
%          res.lambda,res.Cs,res.aucs] = deal(cell(nFolds,1));
%         for iFold = 1:nFolds
%             fprintf('Subject %s: Fold %d of %d.\n',subj,iFold,nFolds)
%             [res.yPred{iFold},...
%                 res.yTest{iFold},...
%                 res.A{iFold},...
%                 res.intercept{iFold},...
%                 res.err{iFold},...
%                 res.lambda{iFold},...
%                 res.Cs{iFold},...
%                 res.aucs{iFold}] = doRegFun(X,Y,folds,iFold,lambda,normType,sessInds);
%         end
%         
%         % calculate metrics
%         perf = mean(vertcat(res.err{:}));
%         res.perf = perf;        
%         yPred = vertcat(res.yPred{:});
%         [~,~,~,res.AUC] = perfcurve(Y,yPred,true);
%         AUC = res.AUC;
%         lambda = mean([res.lambda{:}]);    
%         
%         % create model using full dataset
%         resFull = doRegFullData(X,Y,folds,normType,sessInds);
%     end
    
    subject       = subj;
    params.lambda = lambda;
    if saveOutput
        objLocs = vertcat(events(eventsToUse).objLocs); 
        save(fname,'res','resFull','Y','objLocs','params','perf','tal','AUC');
    end
catch e
    fname = fullfile(saveDir,[subj '_lasso_error.mat']);
    save(fname,'e')
end

function [cBest,tBest,Cs,Ts,auc_pen] = calcPenalty(X,Y,T,folds,sessInds,normType,prctileThresh)
% Calculate the optimal penalty parameter.
%
% Returns: penalty (the optimal C)
%          C (vector C values tried)
%          auc_pen (the auc for each penalty paramter)

Y = double(Y);
Y(Y==0) = -1;
Cs      = logspace(log10(1e-2),log10(1e4),22);
Ts      = 1:length(unique(T));
nCV     = size(folds,1);
auc_pen = NaN(length(Cs),length(Ts));

for thisC = 1:length(Cs)
    thisPen = Cs(thisC);
    
    % set parameters for either L1 or L2 with the current c
    if strcmpi(normType,'L1')
        liblin_param = ['-c ' sprintf('%f',thisPen) ' -s 6 -q'];
    elseif strcmpi(normType,'L2')
        liblin_param = ['-c ' sprintf('%f',thisPen) ' -s 0 -q'];
    end
    
    
    for thisT = 1:length(Ts)       
        tInds = T==Ts(thisT);
        
        % hold test output
        dec_values = [];
        labels = [];
        
        for cv = 1:nCV
            inds = folds(cv,:);
            
            % train data for this cv
            xTrain = X(inds,tInds);
            %         xTrain = [ones(sum(inds),1) xTrain];
            yTrain = Y(inds);
            [xTrain,m,s] = standardize(xTrain,sessInds(inds));
            
            % test data for this cv
            xTest =  X(~inds,tInds);
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
            [pred, acc, dec] = predict(double(yTest),sparse(xTest),model,'-q');
            
            % keep of track of output over all folds
            if model.Label(1) < 0;
                dec = dec * -1;
            end
            dec_values = vertcat(dec_values, dec);
            labels = vertcat(labels, yTest);
        end
        
        % find AUC for this value of C
        [~,~,~,auc_pen(thisC,thisT)] = perfcurve(labels,dec_values,1);
        
    end         
end

% return C and T with first auc above threshold
% [cBest,tBest] = find(auc_pen==max(auc_pen(:)),1,'first');

thresh = prctile(auc_pen(:),prctileThresh);
ind=find(auc_pen>=thresh,1,'first');
[cBest,tBest]=ind2sub(size(auc_pen),ind);
cBest = Cs(cBest);


function [yProbs,yPreds,yTest,A,err,lambda,tBest,Cs,Ts,aucs] = doRegFun(X,Y,T,folds,iFold,lambda,normType,sessInds,events,nestedCvField,prctileThresh,modelEachTime)
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
sessions   = sessInds(trainInds);
yTrainBool = Y(trainInds);
xTrain     = X(trainInds,:);

% if no lambda given, calculate lambda for this fold.
if isempty(lambda)    
    % trainInds  = folds(iFold,setdiff(1:size(folds,2),toRemove));
    % subFolds   = folds(setdiff(1:size(folds,1),iFold),trainInds);    
    [~,subFolds] = createFolds(sessInds(trainInds)',[events(trainInds).(nestedCvField)]);
    [lambda,tBest,Cs,Ts,aucs] = calcPenalty(xTrain,yTrainBool,T,subFolds,sessions,normType,prctileThresh);         
end

% set parameters for either L1 or L2
if strcmpi(normType,'L1')
    liblin_param = ['-c ' sprintf('%f',lambda) ' -s 6 -q'];
elseif strcmpi(normType,'L2')
    liblin_param = ['-c ' sprintf('%f',lambda) ' -s 0 -q'];
end



% make sure the classes are balanced by doing a number of subsampling
% rounds. First see how many observations to remove to balance them
numToRemove = sum(yTrainBool) - sum(~yTrainBool);
maxSamp = 100;
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
    xTrain = X(trainIndsSamp,T==tBest);
    
    % training set y
    yTrain = double(Y(trainIndsSamp));
    yTrain(yTrain==0) = -1;
    
    % standardize training data
    % xTrain = [ones(size(xTrain,1),1) xTrain];
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
       Ws = NaN(100,size(model.w,2)); 
    end
    Ws(nSamp,:) = model.w;
    
    % testing set
    xTest = X(~trainInds,T==tBest);
    yTest = double(Y(~trainInds));
    yTest(yTest==0) = -1;
    
    % standardize testing data
    % xTest = [ones(size(xTest,1),1) xTest];
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


function [nFolds,folds,trials] = createFolds(session,cvField)

[trials,~,trialInds] = unique([session' cvField'],'rows');
nFolds = size(trials,1);
folds = false(nFolds,size(trialInds,1));
for iFold = 1:nFolds
    folds(iFold,:) = trialInds ~= iFold;
end

function res = doRegFullData(X,Y,T,folds,normType,sessInds,prctileThresh)
% train on all data. No test

% find lambda
res = [];
[lambda,bestT,Cs,Ts,aucs] = calcPenalty(X,Y,T,folds,sessInds,normType,prctileThresh); 
res.C = lambda;
res.T = bestT;
res.Cs = Cs;
res.Ts = Ts;
res.aucs = aucs;

% standardize 
[X,m,s] = standardize(X(:,T==bestT),sessInds);

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
sessVec = [events.session];
trialVec = [events.blocknum];

uniqSess = unique(sessVec);
sessPerc = NaN(1,length(sessVec));
for sess = 1:length(uniqSess)   
    sessInd = sessVec == uniqSess(sess);
    sessPerc(sessInd) = [1:sum(sessInd)]./sum(sessInd);            
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














