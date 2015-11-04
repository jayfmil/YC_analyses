function YC1_runMulti_liblinear(subjs,params)
% function YC1_runMulti(subjs,params)
% Inputs:
%
%      subjs - cell array of subject strings (default get_subs('RAM_YC1'))
%     params - params structure (default is returned by multiParams)
%
% Wrapper to YC1_runMulti_subj. By default, will loop over all subjects in
% the YC1 database using the default params returned by multiParams(). This
% performs the lasso regularized logistic regression.
%
% If params.useCorrectedPower is false (default), output is saved in
% params.basePath/OrigPower/<subject>_lasso.mat

% use default params if none given
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

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% see if this was submitted with an open pool. If so, parallel on the level
% of subjects. Otherwise, will loop over subjects one by one.
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    parfor s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        YC1_runMulti_subj_liblinear(subjs{s},params,saveDir);
    end
else    
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        YC1_runMulti_subj_liblinear(subjs{s},params,saveDir);
    end
end









