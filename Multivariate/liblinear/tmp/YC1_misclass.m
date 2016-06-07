function YC1_aucsCorr(subjs,params)
% function YC1_aucsCorr
% use default params if none given
if ~exist('params','var') || isempty(params)
    params = multiParams();
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
        fprintf('Processing %s.\n',params,subjs{s})
        res = load_subj_vars(subjs{s},saveDir);
    end
else    
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        res = load_subj_vars(subjs{s},params,saveDir);
        if ~exist('resAll','var')  && ~isempty(res)
            resAll = res;
        elseif ~isempty(res)
            resAll = merge_structs(resAll,res);
        end
    end
end



keyboard

function res = load_subj_vars(subj,params,saveDir)

res = [];

subjFile = fullfile(saveDir,[subj '_lasso.mat']);
if ~exist(subjFile,'file')
    return
end

% load subject model
subjData = load(subjFile);

% load subject events
events = get_sub_events('RAM_YC1',subj);
events = addErrorField(events);
eventsUsed = events(params.eventFilter(events));

y    = vertcat(subjData.res.yTest{:});
y(y==1) = 2;
y(y==-1) = 1;
pred = vertcat(subjData.res.yPred{:})>0;
pred = double(pred)+1;
probs = vertcat(subjData.res.yProb{:});

[c,order] = confusionmat(y',pred')

% f = @(x) isequal(x(1)>0,x(2)>0);
% mean(cellfun(f,subjData.res.yPred))
% mean(vertcat(subjData.res.err{:}))

keyboard

function elecs = get_regions(tal)

% make sure loctag isn't missing
if ~isfield(tal,'locTag')
    [tal.locTag] = deal('');
end
if sum(cellfun('isempty',{tal.locTag})) == length(tal)
    [tal.locTag] = deal('');
end
missing               = cellfun('isempty',{tal.locTag});
[tal(missing).locTag] = deal('');


% get the electrode indices of brain regions
% locTag based
elecs = [];
elecs.H       = ~cellfun('isempty',regexpi({tal.locTag},['CA1|CA2|CA3|DG|sub']));
elecs.ec      = ~cellfun('isempty',regexpi({tal.locTag},['ec|erc']));
elecs.MTL     = ~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA2|CA3|DG|sub|amy|phc|prc|BA36|erc']));
elecs.MTL     = ~cellfun('isempty',regexpi({tal.locTag},['ec|amy|phc|prc|BA36|erc']));
elecs.ca1     = ~cellfun('isempty',regexpi({tal.locTag},['ca1']));
elecs.ca3     = ~cellfun('isempty',regexpi({tal.locTag},['ca3']));
elecs.dg      = ~cellfun('isempty',regexpi({tal.locTag},['dg']));
elecs.sub     = ~cellfun('isempty',regexpi({tal.locTag},['sub']));
elecs.phc     = ~cellfun('isempty',regexpi({tal.locTag},['phc']));
elecs.prc     = ~cellfun('isempty',regexpi({tal.locTag},['prc']));

% lobe based
elecs.frontal = strcmp({tal.Loc2},'Frontal Lobe');
elecs.occ     = strcmp({tal.Loc2},'Occipital Lobe');
elecs.par     = strcmp({tal.Loc2},'Parietal Lobe');
elecs.temp    = strcmp({tal.Loc2},'Temporal Lobe') & ~elecs.MTL & ~elecs.H;

% new version based on brodmann areas
ba = {tal.Loc5};
%%%% THE SPACE BEFORE THE NUMBER IS IMPORTANT %%%
elecs.aPFC = ~cellfun('isempty',regexpi(ba,[' 10| 11'])) & ~elecs.MTL & ~elecs.H;
elecs.mPFC = ~cellfun('isempty',regexpi(ba,[' 24| 25| 32| 33'])) & ~elecs.MTL & ~elecs.H;
elecs.PFC  = ~cellfun('isempty',regexpi(ba,[' 45| 47| 9| 46'])) & ~elecs.MTL & ~elecs.H;
elecs.TC   = ~cellfun('isempty',regexpi(ba,[' 20| 21| 37'])) & ~elecs.MTL & ~elecs.H;
elecs.PPC  = ~cellfun('isempty',regexpi(ba,[' 7| 40| 39'])) & ~elecs.MTL & ~elecs.H;
elecs.mPC  = ~cellfun('isempty',regexpi(ba,[' 23| 29| 30| 31'])) & ~elecs.MTL & ~elecs.H;
elecs.OC   = ~cellfun('isempty',regexpi(ba,[' 17| 18| 19'])) & ~elecs.MTL & ~elecs.H;


regions = fieldnames(elecs);
elecs.Other = true(1,length(tal));
for r = 1:length(regions)        
    elecs.Other = elecs.Other & ~elecs.(regions{r});    
end


function out = merge_structs(struct1,struct2)
% assums equal fields

fields = fieldnames(struct1);

out = [];
for f = 1:length(fields);
    out.(fields{f}) = [struct1.(fields{f}) struct2.(fields{f})];
end





function events = addErrorField(events)
% add testError field
% add inner field (1 = inner region, 0 = outer region)

testInd = strcmp({events.type},'NAV_TEST');
recEvents = events(testInd);
[events.testError] = deal(NaN);
[events.eucError] = deal(NaN);
[events.testTime] = deal(NaN);
[events.recalled] = deal(NaN);
[events.inner] = deal(NaN);
sessVec = [events.session];
trialVec = [events.blocknum];
for rec = 1:length(recEvents);
    session = recEvents(rec).session;
    trial = recEvents(rec).blocknum;
    err = recEvents(rec).respPerformanceFactor;
    eucErr = recEvents(rec).respDistErr;
    respTime = recEvents(rec).respTravelTime;
    ind = sessVec == session & trialVec == trial;
    [events(ind).testError] = deal(err);
    [events(ind).eucError] = deal(eucErr);
    [events(ind).testTime] = deal(respTime);
    [events(ind).inner] = deal(abs(recEvents(rec).objLocs(1)) < 568/30 && abs(recEvents(rec).objLocs(2)) < 7);
end





























