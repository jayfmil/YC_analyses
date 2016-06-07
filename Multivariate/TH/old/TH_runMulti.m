function TH_runMulti(subjs,params)
% function TH_runMulti(subjs,params)
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

saveDir = params.basePath;
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_TH1');
end

   
for s = 1:length(subjs)
    fprintf('Processing %s.\n',subjs{s})
    TH1_refactor(subjs{s},params,saveDir);
end









