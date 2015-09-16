function YC2_createChanceDists(subjs,params)
% function YC2_applyWeights(subjs,params)
% Inputs:
%
%      subjs - cell array of subject strings (default get_subs('RAM_YC2'))
%     params - params structure (default is returned by multiParams)
%
% Note: in order to run, a subject must also be a RAM_YC1 subject, and
% YC1_runMulti_subj must have been run for that subject.

if ~exist('params','var') || isempty(params)
    params = multiParams();
end

% analysis settings
% -----------------
numIters = 1000;

% save directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
YC1_dir = fullfile(params.basePath,f(params.useCorrectedPower,y));
saveDir = fullfile(YC1_dir,'YC2');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get YC2 subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC2');
end

% see if this was submitted with an open pool. If so, parallel on the level
% of subjects. Otherwise, will loop over subjects one by one.
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    parfor s = 1:length(subjs)
        
        lassoFile  = fullfile(YC1_dir,[subjs{s} '_lasso.mat']);
        errorFile  = fullfile(YC1_dir,[subjs{s} '_error_lasso.mat']);
        chanceFile = fullfile(YC1_dir,[subjs{s} '_chance_perf_dist.mat']);
        if ~exist(lassoFile,'file')
            fprintf('No YC1 lasso data for %s. Skipping.\n',subjs{s})
            continue
        elseif exist(errorFile,'file')
            fprintf('YC1 lasso error file present for %s. Skipping.\n',subjs{s})
            continue
        elseif ~exist(chanceFile,'file')
            fprintf('No YC1 chance file present for %s. Skipping.\n',subjs{s})            
            continue              
        else
            
            % use the same parameters as the real data, but set it to
            % permute the responses
            subjData = load(lassoFile);
            chanceData = load(chanceFile);
            YC1_params = subjData.params;
            YC1_params.powerPath = params.powerPath;
            YC1_params.doPermute = 1;
            YC1_params.saveOutput = 0;
            YC1_params.loadPower = 1;
            
            fname = sprintf('%s_YC2_chance_perf_dist.mat',subjs{s});
            fname = fullfile(saveDir,fname);
            perf_all       = [];
            perfEncAll     = [];
            auc_all        = [];
            aucEncAll      = [];
            aucBest_all    = [];
            aucBestEnc_all = [];
            for i = 1:numIters
                fprintf('Processing %s iteration %d of %d.\n',subjs{s},i,numIters)
                [AUC,AUC_enc,perf,perfEnc,AUC_best,AUC_bestEnc] = YC2_applyWeights(subjs{s},YC1_params,subjData,chanceData,saveDir);
                perf_all          = [perf_all;perf];
                auc_all           = [auc_all;AUC];
                perfEncAll        = [perfEncAll;perfEnc];
                aucEncAll         = [aucEncAll;AUC_enc];
                aucBest_all       = [aucBest_all;AUC_best];
                aucBestEnc_all    = [aucBestEnc_all;AUC_bestEnc];
                if ~isempty(perf_all)
                    parsave(fname,perf_all,auc_all,perfEncAll,aucEncAll,aucBest_all,aucBestEnc_all)
                end
            end
        end
    end
    
elseif isempty(poolobj)
    for s = 1:length(subjs)
        lassoFile  = fullfile(YC1_dir,[subjs{s} '_lasso.mat']);
        errorFile  = fullfile(YC1_dir,[subjs{s} '_error_lasso.mat']);
        chanceFile = fullfile(YC1_dir,[subjs{s} '_chance_perf_dist.mat']);
        if ~exist(lassoFile,'file')
            fprintf('No YC1 lasso data for %s. Skipping.\n',subjs{s})
            continue
        elseif exist(errorFile,'file')
            fprintf('YC1 lasso error file present for %s. Skipping.\n',subjs{s})
            continue
        elseif ~exist(chanceFile,'file')
            fprintf('No YC1 chance file present for %s. Skipping.\n',subjs{s})
            continue
        else
            
            % use the same parameters as the real data, but set it to
            % permute the responses
            subjData = load(lassoFile);
            chanceData = load(chanceFile);
            YC1_params = subjData.params;
            YC1_params.powerPath = params.powerPath;
            YC1_params.doPermute = 1;
            YC1_params.saveOutput = 0;
            YC1_params.loadPower = 1;
            
            fname = sprintf('%s_YC2_chance_perf_dist.mat',subjs{s});
            fname = fullfile(saveDir,fname);
            perf_all   = [];
            perfEncAll = [];
            auc_all    = [];
            aucEncAll  = [];
            for i = 1:numIters
                fprintf('Processing %s iteration %d of %d.\n',subjs{s},i,numIters)
                [AUC,AUC_enc,perf,perfEnc] = YC2_applyWeights(subjs{s},YC1_params,subjData,chanceData,saveDir);
                perf_all   = [perf_all;perf];
                auc_all    = [auc_all;AUC];
                perfEncAll = [perfEncAll;perfEnc];
                aucEncAll  = [aucEncAll;AUC_enc];
                if ~isempty(perf_all)
                    parsave(fname,perf_all,auc_all,perfEncAll,aucEncAll)
                end
            end
        end
    end
end

function parsave(fname,perf_all,auc_all,perfEncAll,aucEncAll,aucBest_all,aucBestEnc_all)
% Because you can't save files directly in parfor loop
save(fname,'perf_all','auc_all','perfEncAll','aucEncAll','aucBest_all','aucBestEnc_all')



