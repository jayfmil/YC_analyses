function YC1_R1059J_wrapper()

basePath = '/scratch/jfm2/YC1/multi/R1006P';
powerPath = '/scratch/jfm2/power8freqs/';

% get list of YC subjects if non given
subjs = {'R1006P'}

% use default params if none given
if ~exist('params','var')
    params = multiParams();
end
params.powerPath = powerPath;
params.crossValStrictness = 1;
params.useKfold = 0;
params.cvField  = 'blocknum';
params.useWatrous = 1;
params.basePath = basePath;
params.region = 'mtl';
params.overwrite = 1;
params.Cs = [];
% loop over subjets

Cs = logspace(log10(1e-2),log10(1e4),22);

timeBins = {};
for i = 1:100
    timeBins{end+1} = [i i];
end
% timeBins = {[65 65]};
% timeBins = {[93 93]}
perfs = NaN(100,22);
aucs  = NaN(100,22);
for s = 1:length(subjs)
    
    % loop over time bin settings
    %     timeBins = {[1 20;21 40;41 60;61 80;81 100],[1 20],[81 100],[1 100]};
    %     timeBins = {[1 1000;1001 2000;2001 3000;3001 4000;4001 5000],...
    %                 [1 1000],...
    %                 [4001 5000],...
    %                 [1 5000]};
    %     fun = @(x) ceil((x+1000)/20);
    %     timeBins = cellfun(fun,timeBins,'uniformoutput',false);
    parfor t = 1:length(timeBins)%65:65%length(timeBins)
        perfs_t = NaN(1,22);
        aucs_t = NaN(1,22);
        t
        % if ~exist('params','var')
        params = multiParams();
        % end
        params.powerPath = powerPath;
        params.crossValStrictness = 1;
        params.useKfold = 0;
        params.cvField  = 'blocknum';
        params.useWatrous = 1;
        params.basePath = basePath;
        params.region = 'mtl';
        params.overwrite = 1;
        params.Cs = [];
        params.timeBins = timeBins{t};
        
        % loop over regularization method
        normTypes = {'L1','L2'};
        for l = 1:1
            params.normType = normTypes{l};
            
            % loop over inner cross val scheme
            nestedCvFields = {'session','blocknum'};
            for cv = 2:2
                params.nestedCvField = nestedCvFields{cv};
                
                % loop over encoding period selection
                encPeriods = {'first','second','combined'};
                for enc = 3:3
                    params.encPeriod = encPeriods{enc};
                    
                    % make save directory
                    pathExt = sprintf('t-%d_l-%d_cv-%d_enc-%d',t,l,cv,enc);
                    saveDir = fullfile(basePath,pathExt);
                    if ~exist(saveDir,'dir')
                        mkdir(saveDir);
                    end
                    
                    % lastly, loop over phase, power, or powre and phase
                    % together (0=power,1=phase,2=both)
                    usePhase = [0,1,2];
                    for p = 1:1
                        params.usePhase = usePhase(p);
                        for lambdas = 1:length(Cs)
                            pathExt = sprintf('Cs_%d',lambdas);
                            saveDirSub = fullfile(saveDir,pathExt);
                            params.Cs = Cs(lambdas);
                            if ~exist(saveDirSub,'dir')
                                mkdir(saveDirSub);
                            end
                            [perfs_t(lambdas),aucs_t(lambdas)] = YC1_watrous_phase(subjs{s},params,saveDir)
%                             ind = sub2ind([100 22],t,lambdas) 
%                             [perfs(ind),aucs(ind)] = YC1_watrous_phase(subjs{s},params,saveDir)
                        end
                    end
                end
            end
        end
        perfs(t,:) = perfs_t
        aucs(t,:) = aucs_t;
    end
end
keyboard






