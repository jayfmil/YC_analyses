function YC1_createAcrossSessFixedParams(subjs)


params = multiParams();
params.cvField = 'session';
params.encPeriod = 'second';
params.timeBins = [251 300];
params.timeBinLabels = {'Wait'};
params.powerPath = '/scratch/jfm2/power8freqs/';
basePath = '/scratch/jfm2/YC1/multi/acrossSess';
params.normType = 'L2';
params.freqBins = [];

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end
subjs = subjs(~strcmp(subjs,'R1062J'));
subjs = subjs(~strcmp(subjs,'R1065J'));

Cs = logspace(log10(1e-6),log10(1e4),22);
aucs = NaN(length(subjs),length(Cs));
aucs_sessAvg = NaN(length(subjs),length(Cs));
for s = 1:length(subjs)
    fprintf('Processing %s.\n',subjs{s})
    for c = 1:length(Cs);
        params.Cs = Cs(c);
        params.basePath = fullfile(basePath,['C',num2str(c)]);
        if ~exist(params.basePath,'dir')
            mkdir(params.basePath)
        end
        clear AUC res
        if exist(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']),'file')
            load(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']));            
        else            
            [perf,AUC,subject,params,res] = YC1_watrous_phase(subjs{s},params,params.basePath);
        end
        if ~isempty(AUC)
            aucs(s,c) = AUC;
            aucs_sessAvg(s,c) = mean([res.fold_auc{:}]);
        end        
    end
end
fname = fullfile(basePath,'aucs_gridSearch_C.mat');
save(fname,'aucs','aucs_sessAvg','subjs','Cs');
keyboard










