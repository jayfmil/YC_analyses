function YC1_runMulti(subjs,anas,ana_func,params,pool)

 
% analysis settings
% -----------------

if ~exist('anas','var') || isempty(anas)
    anas = {};
    ana_funcs = {};
    
    anas{end+1} = 'lassoReg_allEncoding_binary';
    ana_funcs{end+1} = @lassoReg;    
else
    ana_funcs = {str2func(ana_func)};
end

if ~exist('params','var') || isempty(params)
    params = multiParams();
end

for a = 1:length(ana_funcs)
    
    
    % save directory
    f = @(x,y) y{double(x)+1};
    saveDir = fullfile('/data10/scratch/jfm2/YC1/multi',anas{a});
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    
    % get list of YC subjects
    if ~exist('subjs','var') || isempty(subjs)
        subjs = get_subs('RAM_YC1');
    end
    
    ana_func = ana_funcs{a};
    % see if this was submitted with an open pool
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        parfor s = 1:length(subjs)
            fprintf('Processing %s.\n',subjs{s})
            runMulti_subj(subjs{s},ana_func,params,saveDir);
        end
    else
        for s = 1:length(subjs)
            fprintf('Processing %s.\n',subjs{s})
            runMulti_subj(subjs{s},ana_func,params,saveDir);
            
        end
    end
end




function runMulti_subj(subj,ana_func,params,saveDir)

% 
fname = fullfile(saveDir,[subj '_lasso.mat']);
% if exist(fname,'file')
%     return
% end

% load tal structure
tal = getBipolarSubjElecs(subj,1,1,1);
if ~isfield(tal,'locTag') || ~any(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
    fprintf('No MTL electrodes for %s.\n',subj)
    return
end
nElecs = length(tal);

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

% get time and freq settings
freqBins = params.freqBins;
timeBins = params.timeBins;

% load power for all electrodes
% powerData = loadAllPower(tal,subj,events,freqBins,timeBins,config,eventsToUse);
powerData = loadResids_locs(tal,subj,freqBins,timeBins,config,eventsToUse);
powerData = permute(powerData,[3 1 2 4]);
if params.savePower
    powDir = fullfile(saveDir,'power');
    if ~exist(powDir,'dir')
        mkdir(powDir);
    end
    powFile = fullfile(powDir,[subj '_binnedPower.mat']);
    save(powFile,'powerData','params')
end

% determine the cross validation folds. each fold leaves out one trial
[trials,~,trialInds] = unique([session(eventsToUse)' [events(eventsToUse).blocknum]'],'rows');
nFolds = size(trials,1);
folds = false(nFolds,size(trialInds,1));
for iFold = 1:nFolds
    folds(iFold,:) = trialInds ~= iFold;
end

% perform lasso regression on entire data in order to identify best lambda
% value
Y       = [events(eventsToUse).testError]';
YBool   = Y < median(Y);
objLocs = vertcat(events(eventsToUse).objLocs); 

% TO DO: add in regions of env
res = [];
if params.modelEachTime
    for t = 1:size(timeBins,1)
        
        [res(t).yPred,res(t).yTest,res(t).A,res(t).intercept,res(t).err] = deal(cell(nFolds,1));
        for iFold = 1:nFolds
            
            fprintf('Subject %s: Time %d of %d, Fold %d of %d.\n',subj,t,size(timeBins,1),iFold,nFolds)
            X = reshape(squeeze(powerData(:,:,t,:)),size(powerData,1),size(powerData,2)*nElecs);
            [res(t).yPred{iFold},...
                res(t).yTest{iFold},...
                res(t).A{iFold},...
                res(t).intercept{iFold},...
                res(t).err{iFold}] = ana_func(X,Y,YBool,folds(iFold,:));
            
        end
    end
else
    X = reshape(squeeze(powerData),size(powerData,1),size(powerData,2)*size(powerData,3)*nElecs);
    for iFold = 1:nFolds
        fprintf('Subject %s: Fold %d of %d.\n',subj,iFold,nFolds)
        [res.yPred,res.yTest,res.A,res.intercept,res.err] = deal(cell(nFolds,1));
        [res.yPred{iFold},...
            res.yTest{iFold},...
            res.A{iFold},...
            res.intercept{iFold},...
            res.err{iFold}] = ana_func(X,Y,YBool,folds(iFold,:));
    end
    keyboard
end

% % loop over each time bin
% res = [];
% for t = 1:size(timeBins,1)
%     res(t).timeBin = timeBins(t,:);
%     
%     X = reshape(squeeze(powerData(:,:,t,:)),size(powerData,1),size(powerData,2)*nElecs);    
%     fprintf('Calculating best Lambda parameter for timebin %d.\n',t)
%     opt = statset('UseParallel',true);       
%     [B,FitInfo] = lasso(X,Y,'NumLambda',100,'LambdaRatio',.001,'CV',10);
%     lambda = FitInfo.LambdaMinMSE;    
%     
%     res(t).df = FitInfo.DF(FitInfo.IndexMinMSE);
%     res(t).lambda = lambda;
%     res(t).mse  = cell(nFolds,1);
%     res(t).yHat = cell(nFolds,1);
%     res(t).incp = cell(nFolds,1);
%     res(t).beta = cell(nFolds,1);
%     
%     fprintf('Applying..\n')
%     for iFold = 1:nFolds
%         fprintf('Fold %d of %d.\n',iFold,nFolds)
%         [res(t).mse{iFold},res(t).yHat{iFold},res(t).incp{iFold},res(t).beta{iFold}] = ana_func(X,Y,lambda,folds(iFold,:));
%     end       
% end
% save(fname,'res','Y','objLocs','freqBins','timeBins');



% do correct cross validation
% no freq bins
% uniform time bins. try different sizes. Try the whole
% binary on something
% elec x freq x time all at once
% just single sessions?






% set up the weights structure, w/info on class weights
% Weights.MeanCV = mean(WeightMat,2); % mean
% Weights.SECV = nanstd(WeightMat,[],2)/sqrt(size(WeightMat,2)-1);
% Weights.MeanIntercept = mean(W0,2)';

% yHat=vertcat(yHat{:})
% corrThresh = 0:.05:1;
% predThresh = 0:.01:1;
% FA = zeros(length(corrThresh),length(predThresh));
% HR = zeros(length(corrThresh),length(predThresh));
% AUC = zeroes(1,length(corrThresh));
% for t = 1:length(corrThresh)
%     for i = 1:length(predThresh)
%         FA(t,i) = sum((yHat<predThresh(i) & Y > corrThresh(t)))/sum(Y>corrThresh(t));
%         HR(t,i) = sum((yHat<predThresh(i) & Y < corrThresh(t)))/sum(Y<corrThresh(t));
%         AUC(t) = trapz(FA(t,:),HR(t,:));
%     end
% end
% 
% 
% 
% 
% keyboard
function [yPred,yTest,A,intercept,err] = lassoReg(X,Y,YBool,trainInds)

% Under sample larger class?
yTrainBool = YBool(trainInds);

% more recalled than not recalls
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
xHat = mean(xTrain, 2);
xTrain = xTrain - xHat * ones(1,size(xTrain,2));

% training set y
yTrain = Y(trainInds);
yTrain(toRemove) = [];
intercept = mean(yTrain);
yTrain = roundn(yTrain - intercept,-14);
yTrainBool(toRemove) = [];

% compute model
% [A_lasso, stats] = lasso(xTrain', yTrain, 'CV', 5, 'NumLambda', 50);
opt = statset('UseParallel',false);
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    opt = statset('UseParallel',true);
end
[A_lasso, stats] = lassoglm(xTrain',yTrainBool,'binomial','CV', 10, 'NumLambda', 25,'Options',opt);

% testing set
xTest = X(~trainInds,:)';
yTest = Y(~trainInds);
% A = A_lasso(:,stats.IndexMinMSE);
A = A_lasso(:,stats.IndexMinDeviance);
intercept = stats.Intercept(stats.IndexMinDeviance);
B1 = [intercept;A];
yPred = glmval(B1,xTest','logit');
err = mean(round(yPred) == YBool(~trainInds));
% yPred = (xTest - xHat*ones(1,sum(~trainInds)))' * A + intercept;

% mse = mean((yPred-yTest).^2);


%     % create model
%     [B,FitInfo] = lasso(xTrain,yTrain,'lambda',lambda);
%     intercept = FitInfo.Intercept;
%     betas = B;
%
%     % testing set
%     xTest = X(~trainInds,:);
%     yTest = Y(~trainInds);
%
%     % predict
%     yHat = intercept + xTest*betas;
%
%     % compute mse
%     mse = mean((yHat-yTest).^2);



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














