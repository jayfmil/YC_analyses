function YC2_applyWeights(subjs,params)
% function YC2_applyWeights(subjs,params)
% Inputs:
%
%      subjs - cell array of subject strings (default get_subs('RAM_YC1'))
%     params - params structure (default is returned by multiParams)
%
%
% If params.useCorrectedPower is false (default), output is saved in
% params.basePath/OrigPower/<subject>_chance_perf_dist.mat
%
% Note: in order to run

if ~exist('params','var') || isempty(params)
    params = multiParams();
end

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
        if ~exist(lassoFile,'file')
            fprintf('No YC1 lasso data for %s. Skipping.\n',subjs{s})
            continue
        elseif exist(errorFile,'file')
            fprintf('YC1 lasso error file present for %s. Skipping.\n',subjs{s})
            continue
        else
            
            %
            subjData = load(lassoFile);
            YC1_params = subjData.params;
            
            
%             fname = sprintf('%s_chance_perf_dist.mat',subjs{s});
%             fname = fullfile(saveDir,fname);
            
            fprintf('Processing %s.\n',subjs{s})
            YC2_applyWeights(subjs{s},YC1_params,subjData,saveDir);


%             if ~isempty(perf_all)
%                 parsave(fname,perf_all,auc_all)
%             end
        end
    end
end
end
% else
%
%     for s = 1:length(subjs)
%
%         lassoFile  = fullfile(saveDir,[subjs{s} '_lasso.mat']);
%         errorFile  = fullfile(saveDir,[subjs{s} '_error_lasso.mat']);
%         if exist(lassoFile,'file') && ~exist(errorFile,'file')
%
%             % use the same parameters as the real data
%             subjData = load(lassoFile);
%             params = subjData.params;
%             params.doPermute = 1;
%             params.saveOutput = 0;
%             params.loadPower = 1;
%
%             fname = sprintf('%s_chance_perf_dist.mat',subjs{s});
%             fname = fullfile(saveDir,fname);
%             perf_all = [];
%             auc_all = [];
%             for i = 1:numIters
%                 fprintf('Processing %s iteration %d of %d.\n',subjs{s},i,numIters)
%                 [perf,auc] = YC1_runMulti_subj(subjs{s},params,saveDir);
%                 perf_all = [perf_all;perf];
%                 auc_all = [auc_all;auc];
%                 if ~isempty(perf_all)
%                     save(fname,'perf_all','auc_all')
%                 end
%             end
%
%         end
%     end
% end
%
