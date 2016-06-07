function TH_createAcrossTrialFixedC(subjs)


% p = TH_multiParams
% % params.cvField = 'session';
% params.timeBinLabels = {''};
% params.powerPath = '/scratch/jfm2/power8freqs/';
% basePath = '/scratch/jfm2/TH1/multi/acrossTrial';
% params.normType = 'L2';
% params.freqBins = [];

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_TH1');
end


Cs = logspace(log10(1e-6),log10(1e4),22);
aucs = NaN(length(subjs),length(Cs));
for s = 1:length(subjs)
    params = TH_multiParams();
    %     params.cvField = 'session';
    %     params.timeBinLabels = {''};
    params.powerPath = '/scratch/jfm2/power8freqs/';
    params.basePath = '/scratch/jfm2/TH1/multi/acrossTrial';
    params.normType = 'L2';
    params.freqBins = [];
    params.saveOutput = 0;
    
    fprintf('Processing %s.\n',subjs{s})
    for c = 1:length(Cs);
        if c == 1
            params.loadPower=0;
        else
            params.loadPower=1;
        end
        params.Cs = Cs(c);
        
        
        if ~exist(params.basePath,'dir')
            mkdir(params.basePath)
        end
        clear AUC res
        if exist(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']),'file')
            load(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']));
        else
            AUC = TH1_refactor(subjs{s},params,params.basePath);
        end
        if ~isempty(AUC)
            aucs(s,c) = AUC;            
        end
    end
    
end
keyboard
basePath = '/scratch/jfm2/TH1/multi/acrossTrial';
fname = fullfile(basePath,'aucs_gridSearch_C.mat');
save(fname,'aucs','Cs','subjs');
keyboard










