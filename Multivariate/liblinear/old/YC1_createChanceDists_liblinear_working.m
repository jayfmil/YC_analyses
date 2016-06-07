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
numIters = 100;

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
    func = @YC1_runMulti_subj_combinedFeatures_working;
end
    

for s = 1:length(subjs)
    
    lassoFile  = fullfile(saveDir,[subjs{s} '_lasso.mat']);    
    if exist(lassoFile,'file')
        
        % use the same parameters as the real data
        subjData = load(lassoFile);
        params = subjData.params;
        params.doPermute = 1;
        params.saveOutput = 0;
        params.loadPower = 1;        
        
        fname = sprintf('%s_chance_perf_dist_randTrain.mat',subjs{s});
        fname = fullfile(saveDir,fname);
        perf_all = [];
        auc_all = [];        
        for i = 1:numIters
            fprintf('Processing %s iteration %d of %d.\n',subjs{s},i,numIters)
            [perf,auc] = func(subjs{s},params,saveDir);
            perf_all = [perf_all;perf];
            auc_all = [auc_all;auc];            
            if ~isempty(perf_all)
                save(fname,'perf_all','auc_all')
            end
        end
        
    end
end






