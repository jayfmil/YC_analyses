function YC1_searchCVthresh_combinedFeatures(subjs,initialParams,thresholds)
%
%
%

% if isunix && ~ismac
%     p = gcp('nocreate');
%     if isempty(p)
%         open_rhino2_pool(50,'12G');
%     end
% end


% get subjects
% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

if ~exist('initialParams','var') || isempty(initialParams)
    initialParams = multiParams();
end

if ~exist('thresholds','var') || isempty(thresholds)
    thresholds = 75:5:100;
end

parfor s = 1:length(subjs)
    try
        
        % get basic parameters
        params = initialParams;             
%         params.saveOutput = 0;
%         params.cvField = 'session';
        
        % save directory
        f = @(x,y) y{double(x)+1};
        y = {'OrigPower','CorrectedPower'};
        saveDir = fullfile(params.basePath,f(params.useCorrectedPower,y));
        if ~exist(saveDir,'dir')
            mkdir(saveDir);
        end
        
        perfs = NaN(1,length(thresholds));
        aucs = NaN(1,length(thresholds));
        for t = 1:length(thresholds)
            
            if t == 1
                params.savePower = 1;
                params.loadPower = 0;
            else
                params.savePower = 0;
                params.loadPower = 1;
            end 
            
            
            params.auc_prctileThresh = thresholds(t);            
            [perf,auc] = YC1_runMulti_subj_combinedFeatures(subjs{s},params,saveDir);
            if ~isempty(perf)
                perfs(t) = perf;
                aucs(t)  = auc;                         
                fname = fullfile(saveDir,[subjs{s} '_aucs.mat']);
                parsave(fname,aucs,perfs);
            end
        end

    end
end

function parsave(fname,aucs,perfs)
% Because you can't save files directly in parfor loop
save(fname,'aucs','perfs')
