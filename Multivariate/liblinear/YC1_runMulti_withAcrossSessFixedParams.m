function YC1_runMulti_withAcrossSessFixedParams(subjs)


params = multiParams();
params.cvField = 'session';
params.encPeriod = 'second';
params.timeBins = [251 300];
params.timeBinLabels = {'Wait'};
params.powerPath = '/scratch/jfm2/power8freqs/';
params.basePath = '/scratch/jfm2/YC1/multi/acrossSess/gridSearch_C_res';
params.normType = 'L2';
params.freqBins = [];

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end
subjs = subjs(~strcmp(subjs,'R1062J'));
subjs = subjs(~strcmp(subjs,'R1065J'));

% load all subejct results
fname = '/scratch/jfm2/YC1/multi/acrossSess/aucs_gridSearch_C.mat';
stats = load(fname);


for s = 1:length(subjs)
    fprintf('Processing %s.\n',subjs{s})
    
    % choose best C exlcuding res from current subject
    meanAucs = nanmean(stats.aucs(~strcmp(stats.subjs,subjs{s}),:),1);
    C = stats.Cs(find(meanAucs == max(meanAucs),1,'first'));
    params.Cs = C;
    
    if ~exist(params.basePath,'dir')
        mkdir(params.basePath)
    end
    
    [perf,AUC,subject,params,res] = YC1_watrous_phase(subjs{s},params,params.basePath);
end
keyboard







