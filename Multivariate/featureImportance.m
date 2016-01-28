function featureImportance(subjs,params,overwrite)

% if not given, use default params
if ~exist('params','var') || isempty(params)
    params = multiParams();
end

% overwrite exisiting figures?
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

% tex directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
dataDir = fullfile(params.basePath,f(params.useCorrectedPower,y));
saveDir = fullfile(dataDir,'reports');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% figure directory
figDir = fullfile(saveDir,'figs');
if ~exist('figDir','dir')
    mkdir(figDir)
end

featImps = NaN(length(subjs),4);
for s = 1:length(subjs)
    subj = subjs{s};
    fprintf('Creating plots for %s.\n',subj);
    
    % will hold subject specific figure paths
    figs_subj = struct('subj',[],'AUC',[],'recChange',[],'nElecs',[],'bestTime',[]);
    if isempty(params.region)
        params.region = 'all';
    end
    figs_subj.subj   = subj;
    
    % see if files exist for subject. if not, continue
    lassoFile = fullfile(dataDir,[subj '_lasso_pow.mat']);
    if params.usePhase==1
        lassoFile = fullfile(dataDir,[subj '_lasso_phase.mat']);
    elseif params.usePhase==2
        lassoFile = fullfile(dataDir,[subj '_lasso_powphase.mat']);
    end 
    
    if ~exist(lassoFile,'file')
        fprintf('Lasso not found for %s.\n',subj)
        continue
    end
    
    % load subject data
    lassoData  = load(lassoFile);
    events = get_sub_events('RAM_YC1',subj);
    
    % add the test error to the learning trials
    events  = addExtraFields(events);
    session = [events.session];
    
    % filter to events of interest
    eventsToUse = params.eventFilter(events);
    if sum(eventsToUse) < 10
        fprintf('Not enough events for %s.\n',subj)
        return
    end
    encPeriod = [events(eventsToUse).withinItemCount];
    
    %---------------------------------------------------------------------%
    %%% REPLACE ONCE USED POWER/PHASE IS SAVED OUT %%%
    [powerData,sessInds,~,freqLabel] = loadWatrousFeatures(subj,lassoData.params);
    
    
    % make sure the power values correspond to the events I expect
    if sum(eventsToUse) ~= size(powerData,1)
        fprintf('EVENTS MISMATCH FOR %s.\n',subj)
        return
    end
    
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
    %---------------------------------------------------------------------%
    
    %     X         = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs*nItemObs);
    %     T         = reshape(squeeze(timeLabel),size(timeLabel,1),nFreqs*nTimes*nElecs*nItemObs);
    %     trialType = reshape(squeeze(encLabel),size(encLabel,1),nFreqs*nTimes*nElecs*nItemObs);
    
    features = NaN(size(powerData,1),size(powerData,2)*size(powerData,4));
    freqInd  = reshape(squeeze(freqLabel(1,:,1,:)),1,[]);
    
    if size(features,1) ~= length(lassoData.res.A)
        keyboard
    end
    for i = 1:size(features,1)
        
        
        
        timeToUse = lassoData.res.tBest{i};
        encToUse  = lassoData.res.encBest{i};
        
        if encToUse == 3
            feats = squeeze(nanmean(powerData(i,:,timeToUse,:,:),5));
        else
            feats = squeeze(powerData(i,:,timeToUse,:,encToUse));
        end
        
        feats = reshape(feats,1,[]);
        features(i,:) = feats;                        
        
    end
    sessInds = sessInds(encPeriod==1);
    features = standardize(features,sessInds);
    
    
    yProbs = repmat([lassoData.res.yProb{:}]',1,size(features,2));
    bands = unique(freqInd(:));
    featImp = NaN(1,length(bands));
    featImpsAll = NaN(1,size(features,2));
    for f = 1:size(features,2)
        
%         c=cov(features(:,freqInd==b),yProbs(:,freqInd==b));
        c=cov(features(:,f),yProbs(:,f));
        featImpsAll(f) = c(1,2);
    end
    for b = 1:length(bands);
        featImp(b) = nanmean(featImpsAll(freqInd==b));
    end
    featImps(s,:) = featImp
    
end
keyboard




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



function [powerData,sessInds,timeLabel,freqLabel] = loadWatrousFeatures(subj,params)

featureDir = fullfile('/home2/andrew.watrous/Results/YC1_Features',subj,'Features_4bands_all_elecs_11_16_2015');
sessions   = dir(fullfile(featureDir,['Session_*']));


subjDataPow   = cell(1,4);
subjDataPhase = cell(1,4);
sessInds = [];
for s = 1:length(sessions)
    fname = fullfile(featureDir,sessions(s).name,'Subject_features.mat');
    sessData = load(fname);
    mtlElecs = ~cellfun('isempty',{sessData.features.elecs.locTag});
    sessInds = cat(1,sessInds,ones(size(sessData.features.pow.band{1},3),1)*s);
    for band = 1:4
        subjDataPow{band} = cat(3,subjDataPow{band},sessData.features.pow.band{band}(mtlElecs,params.timeBins,:));
        subjDataPhase{band} = cat(3,subjDataPhase{band},sessData.features.phase.band{band}(mtlElecs,params.timeBins,:));
    end
end


dataPow = NaN(length(sessInds),4,size(subjDataPow{1},2),size(subjDataPow{1},1));
dataPhase = NaN(length(sessInds),4*2,size(subjDataPhase{1},2),size(subjDataPhase{1},1));
timeLabelPow = NaN(1,size(dataPow,2),size(dataPow,3),size(dataPow,4));
timeLabelPhase = NaN(1,size(dataPhase,2),size(dataPhase,3),size(dataPhase,4));
freqLabelPow = NaN(1,size(dataPow,2),size(dataPow,3),size(dataPow,4));
freqLabelPhase = NaN(1,size(dataPhase,2),size(dataPhase,3),size(dataPhase,4));

for f = 1:4
    bandPow = subjDataPow{f};
    bandPow = permute(bandPow,[3 2 1]);
    
    bandPhase = subjDataPhase{f};
    bandPhase = permute(bandPhase,[3 2 1]);
    bandSin = sin(bandPhase);
    bandCos = cos(bandPhase);
    
    for t = 1:size(bandPow,2)
        dataPhase(:,f*2-1,t,:) = (bandSin(:,t,:));
        dataPhase(:,f*2,t,:) = (bandCos(:,t,:));
        dataPow(:,f,t,:) = (bandPow(:,t,:));
        timeLabelPow(:,:,t,:) = t;
        timeLabelPhase(:,:,t,:) = t;
        
        freqLabelPow(:,f,:,:) = f;
        freqLabelPhase(:,f*2-1,:,:) = f;
        freqLabelPhase(:,f*2,:,:) = f;
        
    end
    
    
    
end

if params.usePhase==1
    powerData = dataPhase;
    timeLabel = timeLabelPhase;
    freqLabel = freqLabelPhase;
elseif params.usePhase==0
    powerData = dataPow;
    timeLabel = timeLabelPow;
    freqLabel = freqLabelPow;
elseif params.usePhase == 2
    powerData = cat(2,dataPow,dataPhase);
    timeLabel = cat(2,timeLabelPow,timeLabelPhase);
    freqLabel = cat(2,freqLabelPow,freqLabelPhase);
end


tal = sessData.features.elecs(mtlElecs);


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

