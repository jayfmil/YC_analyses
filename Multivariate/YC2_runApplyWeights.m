function YC2_runApplyWeights(subjs,params)
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
            subjData = load(lassoFile);
            YC1_params = subjData.params;
            fprintf('Processing %s.\n',subjs{s})
            YC2_applyWeights(subjs{s},YC1_params,subjData,saveDir);            
        end
    end
elseif isempty(poolobj)
    for s = 1:length(subjs)        
        lassoFile  = fullfile(YC1_dir,[subjs{s} '_lasso.mat']);
        errorFile  = fullfile(YC1_dir,[subjs{s} '_error_lasso.mat']);
        if ~exist(lassoFile,'file')
            fprintf('No YC1 lasso data for %s. Skipping.\n',subjs{s})
            continue
        elseif exist(errorFile,'file')
            fprintf('YC1 lasso error file present for %s. Skipping.\n',subjs{s})
            continue
        else                        
            subjData = load(lassoFile);
            YC1_params = subjData.params;
            YC1_params.powerPath = params.powerPath;
            fprintf('Processing %s.\n',subjs{s})
            YC2_applyWeights(subjs{s},YC1_params,subjData,saveDir);                        
        end
    end
end