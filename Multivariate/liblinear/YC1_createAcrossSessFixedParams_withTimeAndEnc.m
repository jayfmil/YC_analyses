function YC1_createAcrossSessFixedParams_withTimeAndEnc(subjs)


params = multiParams();
params.cvField = 'session';
params.timeBinLabels = {''};
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
Ts = [51 100;101 150;151 200;201 250;251 300;51 300;1 50;301 350;351 400];
encs = {'first','second'};

aucs = NaN(length(subjs),length(Cs),length(Ts),length(encs));
aucs_sessAvg = NaN(length(subjs),length(Cs),length(Ts),length(encs));
parfor s = 1:length(subjs)
    params = multiParams();
    params.cvField = 'session';
    params.timeBinLabels = {''};
    params.powerPath = '/scratch/jfm2/power8freqs/';
    basePath = '/scratch/jfm2/YC1/multi/acrossSess';
    params.normType = 'L2';
    params.freqBins = [];
    fprintf('Processing %s.\n',subjs{s})
    for c = 1:length(Cs);
        params.Cs = Cs(c);
        for t = 1:length(Ts)
            params.timeBins = Ts(t,:);
            for e = 1:length(encs)
                params.encPeriod = encs{e};
                params.basePath = fullfile(basePath,['C',num2str(c),'_T',num2str(t),'_E',num2str(e)]);
                
                if ~exist(params.basePath,'dir')
                    mkdir(params.basePath)
                end
%                 clear AUC res
%                 if exist(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']),'file')
%                     load(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']));
%                 else
                    [perf,AUC,subject,params,res] = YC1_watrous_phase(subjs{s},params,params.basePath);
%                 end
%                 if ~isempty(AUC)
%                     aucs(s,c,t,e) = AUC;
%                     aucs_sessAvg(s,c,t,e) = mean([res.fold_auc{:}]);
%                 end                
            end
        end
    end
end
basePath = '/scratch/jfm2/YC1/multi/acrossSess';
fname = fullfile(basePath,'aucs_gridSearch_C_T_E.mat');
save(fname,'aucs','aucs_sessAvg','subjs','Cs','Ts','encs');
keyboard










