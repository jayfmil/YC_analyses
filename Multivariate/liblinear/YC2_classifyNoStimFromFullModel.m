function res = YC2_classifyNoStimFromFullModel(subjs)

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end
subjs = subjs(~strcmp(subjs,'R1062J'));
subjs = subjs(~strcmp(subjs,'R1065J'));

params = multiParams();
params.cvField = 'session';
params.timeBinLabels = {''};
params.powerPath = '/scratch/jfm2/power8freqs/';
params.normType = 'L2';
params.freqBins = [];

baseDir = '/scratch/jfm2/YC1/multi/acrossSess/fullModels';
saveDir = '/scratch/jfm2/YC1/multi/acrossSess/YC2_noStim';
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

res = [];
parfor s = 1:length(subjs)
       
    fname = fullfile(baseDir,[subjs{s} '_fullModel.mat']);
    if ~exist(fname,'file')
        fprintf('No model for %s. Skipping.\n',subjs{s});
        continue
    end
    
    yc1Data = load(fname);
           
    % run
    try
        resSubj = YC2_classifyNoStimFromFullModel_subj(subjs{s},yc1Data,saveDir);
        res = [res resSubj];
    end    
end
fname = fullfile(saveDir,'all_subj_res.mat');
save(fname,'res');
keyboard

function res = YC2_classifyNoStimFromFullModel_subj(subj,yc1Data,saveDir)
fname = fullfile(saveDir,[subj '_yc2_noStim.mat']);
res = [];

params = yc1Data.params;
% params.encPeriod = 'second'

% load yc1 events to get rec/nonrec thresh
eventsYC1 = get_sub_events('RAM_YC1',subj);
eventsYC1 = addExtraFields(eventsYC1);
yc1thresh = median([eventsYC1(strcmp({eventsYC1.type},'NAV_TEST')).testError]);

% get yc2 events
events = get_sub_events('RAM_YC2',subj);

% add the test error to the learning trials
events  = addExtraFields(events);
session = [events.session];

% filter to non-stim events
evMask = params.eventFilter(events);
eventsToUse = evMask & [events.isStim]~=1;
if sum(eventsToUse) < 48
    fprintf('Not enough events for %s.\n',subj)
    return
end
encPeriod = [events(eventsToUse).withinItemCount];

% load power or phase
if params.useWatrous == 1
    [powerData,sessInds,~,tal] = loadWatrousFeatures(subj,params);
else
    tal = getBipolarSubjElecs(subj,1,1,params.excludeEpiElecs);
    [powerData,sessInds,tal] = loadAllPower(tal,subj,events,params.timeBins,eventsToUse,params);
    powerData = permute(powerData,[3 1 2 4]);
end

if isempty(tal)
    fprintf('no mtl elecs for %s.\n',subj)
    return
end

% make sure the power values correspond to the events I expect
if sum(eventsToUse) ~= size(powerData,1)
    fprintf('EVENTS MISMATCH FOR %s.\n',subj)
    keyboard
    return
end

% filter to just encoding period of intered
if strcmp(params.encPeriod,'first')
    powerData  = powerData(encPeriod==1,:,:,:);
elseif strcmp(params.encPeriod,'second')
    powerData = powerData(encPeriod==2,:,:,:);
end

% size of feature matrix
nElecs = size(powerData,4);
nTimes = size(powerData,3);
nFreqs = size(powerData,2);
nItemObs = size(powerData,5);

% response data
sessInds = sessInds(encPeriod==1);
sessions = session(eventsToUse);
sessions = sessions(encPeriod==1);
if length(sessions) ~= length(sessInds)
    keyboard
end

Y = [events(eventsToUse).testError]';
Y = Y(encPeriod==1);
% Y  = Y < median(Y);
% Y = Y < yc1thresh;

% reshape into obs x features
X   = reshape(squeeze(powerData),size(powerData,1),nFreqs*nTimes*nElecs*nItemObs);

% standardize
[X,m,s] = standardize(X,sessInds);
% Y = double(Y);
% Y(Y==0) = -1;
    
uniqSess = unique(sessInds);
uniqSess = unique(sessions);
for s = 1:length(uniqSess)
    
    thisSess = sessions == uniqSess(s);
    if sum(thisSess) < 24        
        continue
    end
    resSess = [];
    
    % make each session have 50% recalled? Does that even make sense if
    % you're just lloking at nonsitm items?
    ySess = Y(thisSess);
    errs = [events(evMask).testError];
    thresh = median(errs([events(evMask).session]==uniqSess(s)));
%     thresh = median(Y);
    thresh = median(errs);
    
%     ySess = ySess < median(ySess);
    ySess = ySess < thresh;
    ySess = double(ySess);
    ySess(ySess==0) = -1;
    
    % probs = glmval(yc1Data.res.weights',X,'logit','constant','off')';
    [preds,acc,probs] = predict(double(ySess),sparse(X(thisSess,:)),yc1Data.res.model,'-b 1 -q');
    [~,~,~,resSess.auc] = perfcurve(double(ySess),probs(:,1),1);
    resSess.probs = probs;
    
    permAucs = NaN(1,1000);
    for iter = 1:1000 
        randOrder = randperm(length(ySess));
        [preds,acc,probs] = predict(double(ySess(randOrder)),sparse(X(thisSess,:)),yc1Data.res.model,'-b 1 -q');
        [~,~,~,permAucs(iter)] = perfcurve(double(ySess(randOrder)),probs(:,1),1);
    end
    resSess.pval = mean(resSess.auc < permAucs);
    resSess.subj = subj;
    resSess.session = s    
    res = [res resSess];
end
save(fname,'res')


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


function [powerData,sessInds,tal] = loadAllPower(tal,subj,events,timeBins,eventsToUse,params)


nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
[tal] = filterTalByRegion(tal,params.region);
nElecs = length(tal);

powerParams = load(fullfile(params.powerPath,'params_RAM_YC2.mat'));
freqs       = powerParams.params.pow.freqs;

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
    sessions  = unique([events.session]);
    subjPow   = [];
    subjPhase = [];
    sessInds  = [];
    
    for s = 1:length(sessions)
        fname = fullfile(subjPath,'RAM_YC2_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
        sessPow = load(fname);
        
        if e == 1;
            if isempty(params.freqBins)
                nFreqs = size(sessPow.sessOutput.pow,1);
            else
                nFreqs = size(params.freqBins,1);
            end
            nFreqsPhase = size(sessPow.sessOutput.pow,1);
            powerData = NaN(nFreqs,nTimes,nEvents,nElecs);
            phaseData = NaN(nFreqs*2,nTimes,nEvents,nElecs);
        end
        
        % under zscore
        zPow = sessPow.sessOutput.(powField);
        nT = size(sessPow.sessOutput.(powField),2);
        nE = size(sessPow.sessOutput.(powField),3);
        pow = zPow .* repmat(sessPow.sessOutput.(stdField),[1 nT nE]) + repmat(sessPow.sessOutput.(meanField),[1 nT nE]);
        %         pow = zPow;
        
        if nFreqs > 0 && nFreqs ~= length(freqs)
            powTmp = NaN(nFreqs,size(pow,2),size(pow,3));
            for f = 1:nFreqs
                fInds = freqs >= params.freqBins(f,1) & freqs <= params.freqBins(f,2);
                powTmp(f,:,:) = nanmean(pow(fInds,:,:),1);
            end
            pow = powTmp;
        end
        
        subjPow = cat(3,subjPow,pow);
        subjPhase = cat(3,subjPhase,sessPow.sessOutput.phase);
        sessInds = cat(1,sessInds,ones(nE,1)*s);
        
        
    end
    
    if length(eventsToUse) ~= size(subjPow,3)
        fprintf('Number of events does not match size of power matrix for %s!.\n',subj)
        return
    end
    subjPow    = subjPow(:,:,eventsToUse);
    subjPhase  = subjPhase(:,:,eventsToUse);
    sessInds   = sessInds(eventsToUse);
    
    % average times
    for t = 1:nTimes
        tStart = params.timeBins(t,1);
        tEnd = params.timeBins(t,2);
        tInds = tStart:tEnd;
        powerData(:,t,:,e) = nanmean(subjPow(:,tInds,:),2);
        
        % can't specific a dimenstion in circmean, so have to loop?
        if params.usePhase
            for ev = 1:size(phaseData,3);
                for f = 1:nFreqsPhase
                    theta = circmean(subjPhase(f,tInds,ev));
                    phaseData(f*2-1,t,ev,e) = sin(theta);
                    phaseData(f*2,t,ev,e) = cos(theta);
                end
            end
        end
    end
end

if params.usePhase==1
    powerData = phaseData;
elseif params.usePhase == 2
    powerData = cat(1,powerData,phaseData);
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

