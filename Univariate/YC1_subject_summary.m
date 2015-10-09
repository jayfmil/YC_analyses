function YC1_subject_summary(subjs,params)
%function YC1_subject_summary(subjs)
%
%
%


% use default params if none given
if ~exist('params','var') || isempty(params)
    params = univarParams();
end

if size(params.timeBins,1) > 1
    fprintf('Multiple time bins not supported\n')
    return
end

% save directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
region = params.region;
if isempty(region);region = 'all';end
saveDir = fullfile(params.basePath,f(params.useCorrectedPower,y),region);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% see if this was submitted with an open pool. If so, parallel on the level
% of subjects. Otherwise, will loop over subjects one by one.
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    parfor s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        YC1_univarStats(subjs{s},params,saveDir);
    end
else
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        YC1_univarStats(subjs{s},params,saveDir);
    end
end

function YC1_univarStats(subj,params,saveDir)


% load tal structure
tal = getBipolarSubjElecs(subj,params.doBipol,1,params.excludeEpiElecs);
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

% filter to events of interest
eventsToUse = params.eventFilter(events);
er = [events(eventsToUse).testError];
thresh = median(er);
cond1 = [events.testError] < thresh;
cond2 = [events.testError] >= thresh;

% % update the recalled field
% class1 = find([events.testError]<thresh);
% [events(class1).recalled] = deal(1);
% class2 = find([events.testError]>=thresh);
% [events(class2).recalled] = deal(0);


% % conditions of interest
% cond1 = params.ana_func(events, 1);
% cond2 = params.ana_func(events, 0);
% er = [events(cond1|cond2).testError];

if sum(cond1) < 5
    fprintf('Only %d events for %s in cond1 using %s. Skipping subject.\n', sum(cond1),subj,func2str(ana_func))
    return
end

if sum(cond2) < 5
    fprintf('Only %d events for %s in cond2 %s. Skipping subject.\n', sum(cond2),subj,func2str(ana_func))
    return
end

region = params.region;
if isempty(region)
    region = 'all';
end
fprintf('Calculating average power for %d %s elecs.\n',length(tal),region)
powerData = loadAllPower(tal,subj,events,params.freqBins,params.timeBins,powParams,eventsToUse,params);
powerData = permute(powerData,[3 4 1 2]);






% initialize everything
%if isrow(elecs);elecs = elecs';end
% powCond1ByElec = NaN(length(config.distributedParams.freQ),size(elecs,1),'single');
% powCond2ByElec = NaN(length(config.distributedParams.freQ),size(elecs,1),'single');
% powLTA = NaN(size(elecs,1),sum(cond1)+sum(cond2));
% powHTA = NaN(size(elecs,1),sum(cond1)+sum(cond2));
% powG = NaN(size(elecs,1),sum(cond1)+sum(cond2));
% powHFA = NaN(size(elecs,1),sum(cond1)+sum(cond2));
% rLTA =  NaN(size(elecs,1),1);
% pLTA =  NaN(size(elecs,1),1);
% rHTA =  NaN(size(elecs,1),1);
% pHTA =  NaN(size(elecs,1),1);
% rG =  NaN(size(elecs,1),1);
% pG =  NaN(size(elecs,1),1);
% rHFA =  NaN(size(elecs,1),1);
% pHFA =  NaN(size(elecs,1),1);
% statsLTA = struct('tstat',[],'df',[],'sd',[],'p',[],'meanCond1',[],'meanCond2',[]);
% statsHTA = struct('tstat',[],'df',[],'sd',[],'p',[],'meanCond1',[],'meanCond2',[]);
% statsG = struct('tstat',[],'df',[],'sd',[],'p',[],'meanCond1',[],'meanCond2',[]);
% statsHFA = struct('tstat',[],'df',[],'sd',[],'p',[],'meanCond1',[],'meanCond2',[]);
% tagNames = cell(1,size(elecs,1));

% loop over each electrode in region
for e = 1:size(powerData,2)

    for f = 1:size(powerData,3)
        
        pow = powerData(:,e,f)';
        bad = isnan(er) | isnan(pow);
        [r,p] = corr(er(~bad)', pow(~bad)'); 
        [~,p,~,s] = ttest2(pow(cond1),pow(cond2));
        
        
        keyboard
        
%     % corr for low theta
%     test_pow_LTA = nanmean(squeeze(nanmean(pow(fIndLTA,:,cond1|cond2),2)),1);
%     
%     [rLTA(e),pLTA(e)] = corr(er(~bad)', test_pow_LTA(~bad)');
%     powLTA(e,:) = test_pow_LTA;
%     
%     % corr for high theta
%     test_pow_HTA = nanmean(squeeze(nanmean(pow(fIndHTA,:,cond1|cond2),2)),1);
%     bad = isnan(er) | isnan(test_pow_HTA);
%     [rHTA(e),pHTA(e)] = corr(er(~bad)', test_pow_HTA(~bad)');
%     powHTA(e,:) = test_pow_HTA;
%     
%     % corr for gamma
%     test_pow_G = nanmean(squeeze(nanmean(pow(fIndG,:,cond1|cond2),2)),1);
%     bad = isnan(er) | isnan(test_pow_G);
%     [rG(e),pG(e)] = corr(er(~bad)', test_pow_G(~bad)');
%     powG(e,:) = test_pow_G;
%     
%     % corr for HFA
%     test_pow_HFA = nanmean(squeeze(nanmean(pow(fIndHFA,:,cond1|cond2),2)),1);
%     bad = isnan(er) | isnan(test_pow_HFA);
%     [rHFA(e),pHFA(e)] = corr(er(~bad)', test_pow_HFA(~bad)');
%     powHFA(e,:) = test_pow_HFA;
    
    % average across time
    pow = squeeze(nanmean(pow,2));
    
    % t-test cond1 vs cond2 low theta
    [~,p,~,s] = ttest2(nanmean(pow(fIndLTA,cond1)),nanmean(pow(fIndLTA,cond2)));
    statsLTA(e).tstat = s.tstat;
    statsLTA(e).sd = s.sd;
    statsLTA(e).df = s.df;
    statsLTA(e).p = p;
    statsLTA(e).meanCond1 = nanmean(nanmean(pow(fIndLTA,cond1)));
    statsLTA(e).meanCond2 = nanmean(nanmean(pow(fIndLTA,cond2)));
    
    
    % t-test cond1 vs cond2 high theta
    [~,p,~,s] = ttest2(nanmean(pow(fIndHTA,cond1)),nanmean(pow(fIndHTA,cond2)));
    statsHTA(e).tstat = s.tstat;
    statsHTA(e).sd = s.sd;
    statsHTA(e).df = s.df;
    statsHTA(e).p = p;
    statsHTA(e).meanCond1 = nanmean(nanmean(pow(fIndHTA,cond1)));
    statsHTA(e).meanCond2 = nanmean(nanmean(pow(fIndHTA,cond2)));
    
    
    % t-test cond1 vs cond2 gamma
    [~,p,~,s] = ttest2(nanmean(pow(fIndG,cond1)),nanmean(pow(fIndG,cond2)));
    statsG(e).tstat = s.tstat;
    statsG(e).sd = s.sd;
    statsG(e).df = s.df;
    statsG(e).p = p;
    statsG(e).meanCond1 = nanmean(nanmean(pow(fIndLTA,cond1)));
    statsG(e).meanCond2 = nanmean(nanmean(pow(fIndLTA,cond2)));
    
    % t-test cond1 vs cond2 hfa
    [~,p,~,s] = ttest2(nanmean(pow(fIndHFA,cond1)),nanmean(pow(fIndHFA,cond2)));
    statsHFA(e).tstat = s.tstat;
    statsHFA(e).sd = s.sd;
    statsHFA(e).df = s.df;
    statsHFA(e).p = p;
    statsHFA(e).meanCond1 = nanmean(nanmean(pow(fIndHTA,cond1)));
    statsHFA(e).meanCond2 = nanmean(nanmean(pow(fIndHTA,cond2)));
    
    % mean power spect for electrode
    powCond1ByElec(:,e) = nanmean(pow(:,cond1),2);
    powCond2ByElec(:,e) = nanmean(pow(:,cond2),2);
    end
end

% save it to file
save(fname,'powCond1ByElec','powCond2ByElec',...
    'statsLTA','statsHTA','statsG','statsHFA','tagNames',...
    'rLTA','pLTA','rHTA','pHTA','er','powLTA','powHTA',...
    'rG','pG','rHFA','pHFA','powG','powHFA')



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





