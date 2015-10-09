function YC1_subject_summary(subjs,params)
%function YC1_subject_summary(subjs)
 %
 %
 %
 
 
 % use default params if none given
if ~exist('params','var') || isempty(params)
    params = univarParams();
end

% save directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
saveDir = fullfile(params.basePath,f(params.useCorrectedPower,y));
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
thresh = median([events(eventsToUse).testError]);  

% update the recalled field
class1 = find([events.testError]<thresh);
[events(class1).recalled] = deal(1);
class2 = find([events.testError]>=thresh);
[events(class2).recalled] = deal(0);

% LTA freqs
[~,fInd_start] = min(abs(1 - config.distributedParams.freQ));
[~,fInd_end] = min(abs(3 - config.distributedParams.freQ));
fIndLTA = fInd_start:fInd_end;

% HTA freqs
[~,fInd_start] = min(abs(3 - config.distributedParams.freQ));
[~,fInd_end] = min(abs(9 - config.distributedParams.freQ));
fIndHTA = fInd_start:fInd_end;

% GAMMA freqs
[~,fInd_start] = min(abs(40 - config.distributedParams.freQ));
[~,fInd_end] = min(abs(70 - config.distributedParams.freQ));
fIndG = fInd_start:fInd_end;

% HFA freqs
[~,fInd_start] = min(abs(70 - config.distributedParams.freQ));
[~,fInd_end] = min(abs(200 - config.distributedParams.freQ));
fIndHFA = fInd_start:fInd_end;

% conditions of interest
cond1 = ana_func(events, 1);
cond2 = ana_func(events, 0);
er = [events(cond1|cond2).testError];

if sum(cond1) < 5
    fprintf('Only %d events for %s in cond1 using %s. Skipping subject.\n', sum(cond1),subj,func2str(ana_func))
    return
end

if sum(cond2) < 5
    fprintf('Only %d events for %s in cond2 %s. Skipping subject.\n', sum(cond2),subj,func2str(ana_func))
    return
end

% time window of interest
tInds = config.distributedParams.timeBins(:,1) > 1 & ...
        config.distributedParams.timeBins(:,2) <= 4000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TOM: Here is is set to only do hippocampus %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
% load power for a region
for roi = {'hipp'};%,'ec','mtl','frontal','parietal','temporal','occipital','limbic'}
        
    % filter electrodes by current region
    switch roi{1}
        case 'hipp'
            elecs = hipp_elecs;
            if isempty(elecs)
                fprintf('No hipp electrodes for %s.\n',subj)
                continue
            end
            fprintf('Calculating average power for %d hipp elecs.\n',length(elecs))
            fname = fullfile(saveDir,[subj,'_hipp_pow']);
            saveHippDist = true;
        case 'ec'
            elecs = ec_elecs;
            if isempty(elecs)
                fprintf('No EC electrodes for %s.\n',subj)
                continue
            end            
            fprintf('Calculating average power for %d EC elecs.\n',length(elecs))            
            fname = fullfile(saveDir,[subj,'_ec_pow']);
        case 'mtl'
            elecs = mtl_elecs;
            if isempty(elecs)
                fprintf('No MTL electrodes for %s.\n',subj)
                continue
            end            
            fprintf('Calculating average power for %d MTL elecs.\n',length(elecs))            
            fname = fullfile(saveDir,[subj,'_mtl_pow']);
        case 'frontal'
            elecs = frontal_elecs;
            if isempty(elecs)
                fprintf('No frontal surface electrodes for %s.\n',subj)
                continue
            end            
            fprintf('Calculating average power for %d frontal elecs.\n',length(elecs))            
            fname = fullfile(saveDir,[subj,'_frontal_pow']);
        case 'temporal'
            elecs = temporal_elecs;
            if isempty(elecs)
                fprintf('No temporal surface electrodes for %s.\n',subj)
                continue
            end            
            fprintf('Calculating average power for %d temporal elecs.\n',length(elecs))            
            fname = fullfile(saveDir,[subj,'_temporal_pow']);
        case 'occipital'
            elecs = occipital_elecs;
            if isempty(elecs)
                fprintf('No surface occipital electrodes for %s.\n',subj)
                continue
            end            
            fprintf('Calculating average power for %d occipital elecs.\n',length(elecs))            
            fname = fullfile(saveDir,[subj,'_occipital_pow']);
        case 'parietal'
            elecs = parietal_elecs;
            if isempty(elecs)
                fprintf('No surface parietal electrodes for %s.\n',subj)
                continue
            end            
            fprintf('Calculating average power for %d parietal elecs.\n',length(elecs))            
            fname = fullfile(saveDir,[subj,'_parietal_pow']);          
        case 'limbic'
            elecs = limbic_elecs;
            if isempty(elecs)
                fprintf('No surface limbic electrodes for %s.\n',subj)
                continue
            end            
            fprintf('Calculating average power for %d limbic elecs.\n',length(elecs))            
            fname = fullfile(saveDir,[subj,'_limbic_pow']);                    
        otherwise
            disp('UNKNOWN ROI')
            fname = [];
            continue
    end                   
    
    % initialize everything
    %if isrow(elecs);elecs = elecs';end
    powCond1ByElec = NaN(length(config.distributedParams.freQ),size(elecs,1),'single');
    powCond2ByElec = NaN(length(config.distributedParams.freQ),size(elecs,1),'single');    
    powLTA = NaN(size(elecs,1),sum(cond1)+sum(cond2));    
    powHTA = NaN(size(elecs,1),sum(cond1)+sum(cond2));    
    powG = NaN(size(elecs,1),sum(cond1)+sum(cond2));    
    powHFA = NaN(size(elecs,1),sum(cond1)+sum(cond2));    
    rLTA =  NaN(size(elecs,1),1);
    pLTA =  NaN(size(elecs,1),1);
    rHTA =  NaN(size(elecs,1),1);
    pHTA =  NaN(size(elecs,1),1);
    rG =  NaN(size(elecs,1),1);
    pG =  NaN(size(elecs,1),1);
    rHFA =  NaN(size(elecs,1),1);
    pHFA =  NaN(size(elecs,1),1);
    statsLTA = struct('tstat',[],'df',[],'sd',[],'p',[],'meanCond1',[],'meanCond2',[]);
    statsHTA = struct('tstat',[],'df',[],'sd',[],'p',[],'meanCond1',[],'meanCond2',[]);
    statsG = struct('tstat',[],'df',[],'sd',[],'p',[],'meanCond1',[],'meanCond2',[]);
    statsHFA = struct('tstat',[],'df',[],'sd',[],'p',[],'meanCond1',[],'meanCond2',[]);
    tagNames = cell(1,size(elecs,1));
    
    % loop over each electrode in region
    for e = 1:size(elecs,1)
        elecNum = elecs(e,:);
        tagNames{e} = tal(ismember(vertcat(tal.channel),elecNum,'rows')).tagName;

        % load power for all sessions. Power should aleady have been
        % created or else error
        if ~useResids
            pow  = loadPow_local(subj,elecNum,config,events);
        else
            if e == 1     
                eventsToUse = cond1|cond2;
                cond1 = ana_func(events(eventsToUse), 1);
                cond2 = ana_func(events(eventsToUse), 0);                            
            end
            pow = loadResids_locs(subj,elecNum,cond1|cond2);
        end
              
        % use only time bins of interest
        pow(:,~tInds,:) = NaN;
                
        % corr for low theta
        test_pow_LTA = nanmean(squeeze(nanmean(pow(fIndLTA,:,cond1|cond2),2)),1);        
        bad = isnan(er) | isnan(test_pow_LTA);
        [rLTA(e),pLTA(e)] = corr(er(~bad)', test_pow_LTA(~bad)');
        powLTA(e,:) = test_pow_LTA;

        % corr for high theta
        test_pow_HTA = nanmean(squeeze(nanmean(pow(fIndHTA,:,cond1|cond2),2)),1);
        bad = isnan(er) | isnan(test_pow_HTA);
        [rHTA(e),pHTA(e)] = corr(er(~bad)', test_pow_HTA(~bad)');
        powHTA(e,:) = test_pow_HTA;

        % corr for gamma
        test_pow_G = nanmean(squeeze(nanmean(pow(fIndG,:,cond1|cond2),2)),1);
        bad = isnan(er) | isnan(test_pow_G);
        [rG(e),pG(e)] = corr(er(~bad)', test_pow_G(~bad)');
        powG(e,:) = test_pow_G;

        % corr for HFA
        test_pow_HFA = nanmean(squeeze(nanmean(pow(fIndHFA,:,cond1|cond2),2)),1);
        bad = isnan(er) | isnan(test_pow_HFA);
        [rHFA(e),pHFA(e)] = corr(er(~bad)', test_pow_HFA(~bad)');
        powHFA(e,:) = test_pow_HFA;

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
    
    % save it to file
        save(fname,'powCond1ByElec','powCond2ByElec',...            
            'statsLTA','statsHTA','statsG','statsHFA','tagNames',...
             'rLTA','pLTA','rHTA','pHTA','er','powLTA','powHTA',...
             'rG','pG','rHFA','pHFA','powG','powHFA')
    
end

function pow = loadPow_local(subj,elecNum,config,events)
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
        subjMean = nanmean(squeeze(nanmean(subjPow(:,:,inds), ...
            2)),2);
        subjMean = repmat(subjMean,[1 size(subjPow,2), size(subjPow,3)]);
        
        subjStd = nanstd(squeeze(nanmean(subjPow(:,:,inds), ...
            2)),[],2);
        subjStd = repmat(subjStd,[1 size(subjPow,2), size(subjPow,3)]);
        zpow(:,:,inds) = (subjPow(:,:,inds) - subjMean).*subjStd;
        
    end
    subjPow = zpow;
end

% if this happens, the events and power do not correspond
if size(subjPow,3) ~= length(events)
    keyboard
end

% replace time periods outside of each event with nans
pow = subjPow;

function [pow] = loadResids_locs(subj,elecNum,eventsToUse)

basePath  = '/data10/scratch/jfm2/YC1/multi/power/regress/';
subjPath  = fullfile(basePath,subj);
fname     = sprintf('%s_elec_%d-%d_residuals.mat',subj,elecNum(1),elecNum(2));

if ~exist(fullfile(subjPath,fname),'file')
    error('Residuals file %s not found.\n',fname)
else
    elecData = load(fullfile(subjPath,fname));
    pow = permute(elecData.resid,[3 2 1]);
    if size(elecData.resid,1) ~= sum(eventsToUse)
        keyboard
    end
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





