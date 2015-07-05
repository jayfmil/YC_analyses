function YC1_subject_summary(subjs,anas,ana_func,pool)

 
% analysis settings
% -----------------

if ~exist('anas','var') || isempty(anas)

    anas = {};
    ana_funcs = {};
    
    anas{end+1} = 'correct_incorrect_resids';
    ana_funcs{end+1} = @correctFilter;
    
else
    ana_funcs = {str2func(ana_func)};
end


for a = 1:length(ana_funcs)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TOM: Path is set here %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save directory
    f = @(x,y) y{double(x)+1};
    saveDir = fullfile('/data10/scratch/jfm2/YC1/group',anas{a});
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TOM: I don't know if switching this to bipolar will work %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % do bipolar
    bipol = 1;
    
    % use residuals from regression
    useResids = 1;
    
    % get list of YC subjects
    if ~exist('subjs','var') || isempty(subjs)
        subjs = get_subs('RAM_YC1');
    end

    ana_func = ana_funcs{a};

    if exist('pool','var')
        matlabpool(pool,length(subjs)+1)        
        tic
        parfor s = 1:length(subjs)
            fprintf('Processing %s.\n',subjs{s})
            save_pow_in_region(subjs{s},bipol,useResids,ana_func,saveDir);
        end
        toc
        matlabpool close
        toc
    else
        for s = 1:length(subjs)
            fprintf('Processing %s.\n',subjs{s})
            save_pow_in_region(subjs{s},bipol,useResids,ana_func,saveDir);

        end
    end
end

function save_pow_in_region(subj,bipol,useResids,ana_func,saveDir)

% load tal structure
tal = getBipolarSubjElecs(subj,bipol,1);

[hipp_elecs, ec_elecs, mtl_elecs] = deal([]);

% load subject localizations. This info should also be in the tal struct
% now
locs = getLocalization(subj);
if ~isfield(locs,'tag')
  fprintf('localization missing for %s\n',subj)
  return
end

if bipol
  pairs = cellfun('length',{locs.channel})==2;
  locs = locs(pairs);
end

% filter for hipp
hipp_tal = locs(~cellfun('isempty',regexpi({locs.tag},['hipp|CA1|CA3|DG|sub'])));
hipp_elecs = vertcat(hipp_tal.channel);

% filter for ec
ec_tal = locs(~cellfun('isempty',regexpi({locs.tag},['ec'])));
ec_elecs = vertcat(ec_tal.channel);

% filter for all mtl
mtl_tal = locs(~cellfun('isempty',regexpi({locs.tag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])));
mtl_elecs = vertcat(mtl_tal.channel);

% also do lobes
% surface_elecs = tal(~strcmp({tal.eType},'D'));
other_elecs = tal(~ismember(vertcat(tal.channel),[hipp_elecs; ec_elecs; mtl_elecs],'rows'));
lobes = {other_elecs.Loc2};
temporal_elecs = vertcat(other_elecs(strcmp(lobes, 'Temporal Lobe')).channel);
occipital_elecs = vertcat(other_elecs(strcmp(lobes, 'Occipital Lobe')).channel);
frontal_elecs = vertcat(other_elecs(strcmp(lobes, 'Frontal Lobe')).channel);
parietal_elecs = vertcat(other_elecs(strcmp(lobes, 'Parietal Lobe')).channel);
limbic_elecs = vertcat(other_elecs(strcmp(lobes, 'Limbic Lobe')).channel);

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
testInd = strcmp({events.type},'NAV_TEST');
recEvents = events(testInd);
[events.testError] = deal(NaN);
[events.recalled] = deal(NaN);
sessVec = [events.session];
trialVec = [events.blocknum];
for rec = 1:length(recEvents);
  session = recEvents(rec).session;  
  trial = recEvents(rec).blocknum;
  err = recEvents(rec).respPerformanceFactor;
  ind = sessVec == session & trialVec == trial;
  [events(ind).testError] = deal(err);
end

% filter between good and bad
thresh = median([recEvents.respPerformanceFactor]);

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
            [pow,cond1,cond2] = loadResids_locs(subj,elecNum,events,cond1|cond2);
        end
        keyboard
       
        % use only time bins of interest
        pow(:,~tInds,:) = NaN;
                
        % corr for low theta
        test_pow_LTA = nanmean(squeeze(nanmean(pow(fIndLTA,:,cond1|cond2),2)),1);        
        bad = isnan(err) | isnan(test_pow_LTA);
        [rLTA(e),pLTA(e)] = corr(er(~bad)', test_pow_LTA(~bad)');
        powLTA(e,:) = test_pow_LTA;

        % corr for high theta
        test_pow_HTA = nanmean(squeeze(nanmean(pow(fIndHTA,:,cond1|cond2),2)),1);
        bad = isnan(err) | isnan(test_pow_HTA);
        [rHTA(e),pHTA(e)] = corr(er(~bad)', test_pow_HTA(~bad)');
        powHTA(e,:) = test_pow_HTA;

        % corr for gamma
        test_pow_G = nanmean(squeeze(nanmean(pow(fIndG,:,cond1|cond2),2)),1);
        bad = isnan(err) | isnan(test_pow_G);
        [rG(e),pG(e)] = corr(er(~bad)', test_pow_G(~bad)');
        powG(e,:) = test_pow_G;

        % corr for HFA
        test_pow_HFA = nanmean(squeeze(nanmean(pow(fIndHFA,:,cond1|cond2),2)),1);
        bad = isnan(err) | isnan(test_pow_HFA);
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

function [pow,cond1,cond2] = loadResids_locs(subj,elecNum,events,eventsToUse)

basePath  = '/data10/scratch/jfm2/YC1/multi/power/regress/';
subjPath  = fullfile(basePath,subj);
fname     = sprintf('%s_elec_%d-%d_residuals.mat',subj,elecNum(1),elecNum(2));

if ~exist(fullfile(subjPath,fname),'file')
    error('Residuals file %s not found.\n',fname)
else
    elecData = load(fullfile(subjPath,fname));
    if size(elecData.resid,1) ~= sum(eventsToUse)
        keyboard
    end
    pow = permute(elecData.resid,[3 2 1]);    
    cond1 = ana_func(events(eventsToUse), 1);
    cond2 = ana_func(events(eventsToUse), 0);
end


function doNothing(varargin)
% GIVE ERROR. RAM_loadPow requires a power creation function. If this
% code is executed, it is because the power for the subject/session hasn't
% been computed yet. I don't want to compute power here.

error('POWER NOT COMPUTED YET')


function eventsMask = correctFilter(events, isCorrect)     

encInd = strcmp({events.type},'NAV_LEARN') | strcmp({events.type}, ...
                                                  'NAV_LEARN_HARD');
if ~isCorrect
  eventsMask = [events.recalled]==0 & encInd;
else
  eventsMask = [events.recalled]==1 & encInd;
end




