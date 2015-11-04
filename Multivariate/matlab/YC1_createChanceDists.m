function YC1_createChanceDists(subjs,params,justPermuteResponses)
% function YC1_createChanceDists(subjs,params)
% Inputs:
%
%                  subjs - cell array of subject strings (default get_subs('RAM_YC1'))
%                 params - params structure (default is returned by multiParams)
%   justPermuteResponses - If false, will redo classifcation on permuted
%                          response labels. If true, will not redo
%                          classifation, but rather will shuffle the
%                          responses and recalculate AUC and performance
%                          measures. These test slighyl different
%                          hypothesis.. When false, you are asking, "can we
%                          create a model where brain signals fit the
%                          actual behavior better than we can create a
%                          model where brain signals are fit to permuted
%                          responses?" When true, you are asking, "do the
%                          predications of the model fit the real behavior
%                          better than they fit permuted responses?"
%
% Wrapper to YC1_runMulti_subj, but unlike YC1_runMulti, the parameter to
% permute the predicted variable is set to true and this is repeated 1000
% times. The goal is to create a null distribution of classifier accuracy
% to which we can compare the true accuracy.
%
% If params.useCorrectedPower is false (default), output is saved in
% params.basePath/OrigPower/<subject>_chance_perf_dist.mat
%
% Note: in order to run either YC1_makeSubjectReports and
% YC1_weightsByRegion, the chance distributions must exist.

% analysis settings
% -----------------
numIters = 1000;

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

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

if ~exist('justPermuteResponses','var') || isempty(justPermuteResponses)
    justPermuteResponses = 0;
end

if justPermuteResponses
    func = @YC1_permuteOutput;
else
    func = @YC1_runMulti_subj;
end
    

% see if this was submitted with an open pool. If so, parallel on the level
% of subjects. Otherwise, will loop over subjects one by one.
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    parfor s = 1:length(subjs)
        
        lassoFile  = fullfile(saveDir,[subjs{s} '_lasso.mat']);
        errorFile  = fullfile(saveDir,[subjs{s} '_error_lasso.mat']); ...            
        chanceFile = fullfile(saveDir,[subjs{s} '_chance_perf_dist.mat']);
        if exist(lassoFile,'file') && ~exist(errorFile,'file') && ~exist(chanceFile,'file')
            
            % use the same parameters as the real data, but set it to
            % permute the responses
            subjData = load(lassoFile);
            params = subjData.params;
            params.doPermute = 1;
            params.saveOutput = 0;            
            params.loadPower = 1;
%             params.powerPath = '/data10/scratch/jfm2/RAM/biomarker/power/';

            fname = sprintf('%s_chance_perf_dist.mat',subjs{s});
            fname = fullfile(saveDir,fname);
            perf_all = [];
            auc_all  = [];
            r_all   = [];
            for i = 1:numIters
                fprintf('Processing %s iteration %d of %d.\n',subjs{s},i,numIters)
                [perf,auc,r] = func(subjs{s},params,saveDir);
                perf_all = [perf_all;perf];
                auc_all = [auc_all;auc];
                r_all = [r_all;r];
                if ~isempty(perf_all)
                    parsave(fname,perf_all,auc_all,r_all)
                end
            end
        end
    end
else
    
    for s = 1:length(subjs)
        
        lassoFile  = fullfile(saveDir,[subjs{s} '_lasso.mat']);
        errorFile  = fullfile(saveDir,[subjs{s} '_error_lasso.mat']);
        if exist(lassoFile,'file') && ~exist(errorFile,'file')
            
            % use the same parameters as the real data
            subjData = load(lassoFile);
            params = subjData.params;
            params.doPermute = 1;     
            params.saveOutput = 0;
            params.loadPower = 1;
%             params.powerPath = '/data10/scratch/jfm2/RAM/biomarker/power/';

            fname = sprintf('%s_chance_perf_dist.mat',subjs{s});
            fname = fullfile(saveDir,fname);
            perf_all = [];
            auc_all = [];
            r_all   = [];
            for i = 1:numIters
                fprintf('Processing %s iteration %d of %d.\n',subjs{s},i,numIters)
                [perf,auc,r] = func(subjs{s},params,saveDir);
                perf_all = [perf_all;perf];
                auc_all = [auc_all;auc];
                r_all = [r_all;r];
                if ~isempty(perf_all)
                    save(fname,'perf_all','auc_all','r_all')
                end
            end
            
        end
    end
end


function parsave(fname,perf_all,auc_all,r_all)
% Because you can't save files directly in parfor loop
save(fname,'perf_all','auc_all','r_all')







