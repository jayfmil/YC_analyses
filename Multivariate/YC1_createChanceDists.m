function YC1_createChanceDists(subjs,params)


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


% see if this was submitted with an open pool
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    parfor s = 1:length(subjs)
        
        lassoFile  = fullfile(saveDir,[subjs{s} '_lasso.mat']);
        errorFile  = fullfile(saveDir,[subjs{s} '_error_lasso.mat']);        
        if exist(lassoFile,'file') && ~exist(errorFile,'file')
            
            % use the same parameters as the real data
            subjData = load(lassoFile);
            params = subjData.params;
            params.doPermute = 1;
            params.saveOutput = 0;            
            
            fname = sprintf('%s_chance_perf_dist.mat',subjs{s});
            fname = fullfile(saveDir,fname);
            perf_all = [];
            for i = 1:numIters
                fprintf('Processing %s iteration %d of %d.\n',subjs{s},i,numIters)
                perf = YC1_runMulti_subj(subjs{s},params,saveDir);
                perf_all = [perf_all;perf]
                if ~isempty(perf_all)
                    parsave(fname,perf_all)
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
            
            fname = sprintf('%s_chance_perf_dist.mat',subjs{s});
            fname = fullfile(saveDir,fname);
            perf_all = [];
            for i = 1:numIters
                fprintf('Processing %s iteration %d of %d.\n',subjs{s},i,numIters)
                perf = YC1_runMulti_subj(subjs{s},params,saveDir);
                perf_all = [perf_all;perf]
                if ~isempty(perf_all)
                    save(fname,'perf_all')
                end
            end
            
        end
    end
end


function parsave(fname,perf_all)
% Because you can't save files directly in parfor loop
save(fname,'perf_all')







