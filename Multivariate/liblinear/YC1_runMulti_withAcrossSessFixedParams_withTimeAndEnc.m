function YC1_runMulti_withAcrossSessFixedParams_withTimeAndEnc(subjs)


params = multiParams();
params.cvField = 'session';
% params.encPeriod = 'second';
% params.timeBins = [251 300];
% params.timeBinLabels = {'Wait'};
params.powerPath = '/scratch/jfm2/power8freqs/';
params.basePath = '/scratch/jfm2/YC1/multi/acrossSess/gridSearch_C_T_E_res';
params.normType = 'L2';
params.freqBins = [];
params.overwrite = 1;

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end
subjs = subjs(~strcmp(subjs,'R1062J'));
subjs = subjs(~strcmp(subjs,'R1065J'));

% load all subejct results
fname = '/scratch/jfm2/YC1/multi/acrossSess/aucs_gridSearch_C_T_E.mat';
stats = load(fname);


for s = 1:length(subjs)
    fprintf('Processing %s.\n',subjs{s})
    
    
    % choose best C, time, and enc period excluding results from current subject
    meanAucs = squeeze(nanmean(stats.aucs(~strcmp(stats.subjs,subjs{s}),:,:,:),1));
    [cInd,tInd,eInd] = ind2sub(size(meanAucs),find(meanAucs == max(meanAucs(:)),1,'first'));
    params.Cs        = stats.Cs(cInd);
    params.timeBins  = stats.Ts(tInd,:);
    params.encPeriod = stats.encs{eInd};
    
    
    if ~exist(params.basePath,'dir')
        mkdir(params.basePath)
    end
    
    [perf,AUC,subject,params,res] = YC1_watrous_phase(subjs{s},params,params.basePath);
    keyboard
end








