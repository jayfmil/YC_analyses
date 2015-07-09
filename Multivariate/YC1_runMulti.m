function YC1_runMulti(subjs,params)


% analysis settings
% -----------------

if ~exist('params','var') || isempty(params)
    params = multiParams();
end

keyboard

% save directory
f = @(x,y) y{double(x)+1};
saveDir = fullfile('/data10/scratch/jfm2/YC1/multi',ana_name);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end


% see if this was submitted with an open pool
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    parfor s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        YC1_runMulti_subj(subjs{s},params,saveDir);
    end
else
    perf_all = [];
    subj_all = {};
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        [perf,subj] = YC1_runMulti_subj(subjs{s},params,saveDir);
        perf_all = [perf_all;perf]
        subj_all = vertcat(subj_all,subj)
    end
end


keyboard








