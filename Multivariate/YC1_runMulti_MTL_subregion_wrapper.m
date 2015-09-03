function YC1_runMulti_MTL_subregion_wrapper
%
% Run classification seperately for every MTL subregion.
%

if isunix && ~ismac
    p = gcp('nocreate');
    if isempty(p)
        open_rhino2_pool(35,'12G');
    end
end
basePath = '/data10/scratch/jfm2/YC1/multi/mtl_subregions'
timeBinLabels = {'Spin', 'Drive', 'Wait'};
timeBinLabels = {'Pre'  'Spin'  'Drive1'  'Drive2'  'Drive3'  'Wait'  'Post'};
% get subjects
subjs = get_subs('RAM_YC1');
% subjs = subjs(~strcmp(subjs,'R1001P'));

% loop over mtl regions
regions = {'CA1','CA3','DG','EC','PHC','PRC','Sub'};
combs   = nchoosek(1:length(regions),2);

parfor s = 1:length(subjs)
    
    % get basic parameters
    params = multiParams();
    params.modelEachTime = 0;    
%     timeBins = [1     1000;...
%         1001  4000;...        
%         4001  5000];   
    
    timeBins = [-999 0;...
    1 1000;...
    1001 2000;...
    2001 3000;...
    3001 4000;...
    4001 5000;...
    5001 6000];
        
    params.timeBins = timeBins;
    params.timeBinLabels = timeBinLabels;
    params.savePower = 0;
    params.nCV = 10;
    params.excludeEpiElecs = 1;
    
    % first do all MTL
    params.basePath = fullfile(basePath,'MTL');
    params.region = 'MTL';
    YC1_runMulti({subjs{s}},params);
    
    % then loop over each region
    for r = 1:length(regions)
        params.basePath = fullfile(basePath,regions{r});
        params.region = regions{r};
        YC1_runMulti({subjs{s}},params);    
    end
    
    % then loop over all combinations of regions
    for c = 1:size(combs,1)                       
        
        % stupid try here because getBipolarSubjElecs errors if the tal
        % file isn't there
        try
            % make sure subject has electrodes in both regions
            tal  = getBipolarSubjElecs(subjs{s},1,1,params.excludeEpiElecs);
            tal1 = filterTalByRegion(tal,regions{combs(c,1)});
            tal2 = filterTalByRegion(tal,regions{combs(c,2)});
            if isempty(tal1) || isempty(tal2)
                fprintf('%s does not have electrodes in both %s and %s.\n',subjs{s},regions{combs(c,1)},regions{combs(c,2)});
                continue
            end
            
            regionPath = sprintf('%s-%s',regions{combs(c,1)},regions{combs(c,2)});
            regionStr  = sprintf('%s|%s',regions{combs(c,1)},regions{combs(c,2)});
            params.basePath = fullfile(basePath,regionPath);
            params.region = regionStr;
            YC1_runMulti({subjs{s}},params);
        end
    end
end