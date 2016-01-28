function YC1_runMulti_gridWrapper_fixed(subjs,params)

dataPath = '/scratch/jfm2/YC1/multi/pow8_gridWrapper';
basePath = '/scratch/jfm2/YC1/multi/pow8_gridWrapper_fixed_block';
powerPath = '/scratch/jfm2/power8freqs/';

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% use default params if none given
% if ~exist('params','var') || isempty(params)
%     params = multiParams();
% end

% loop over subjets
parfor s = 1:length(subjs)
    
    params = multiParams();
    params.powerPath = powerPath;
    params.crossValStrictness = 1;
    params.useKfold = 0;
    params.cvField  = 'session';
    params.useWatrous = 0;
    params.basePath = basePath;
    
    % loop over time bin settings
    %     timeBins = {[1 20;21 40;41 60;61 80;81 100],[1 20],[81 100],[1 100]};
    timeBins = {[1 1000;1001 2000;2001 3000;3001 4000;4001 5000],...
                [1 1000],...
                [4001 5000],...
                [1 5000]};
    fun = @(x) ceil((x+1000)/20);
    timeBins = cellfun(fun,timeBins,'uniformoutput',false);
    for t = 1:length(timeBins)
        params.timeBins = timeBins{t};
        
        % loop over regularization method
        normTypes = {'L1','L2'};
        for l = 1:2
            params.normType = normTypes{l};
            
            % loop over inner cross val scheme
            nestedCvFields = {'session','blocknum'};
            for cv = 2:2
                params.nestedCvField = nestedCvFields{cv};
                
                % loop over encoding period selection
                encPeriods = {'first','second','combined'};
                for enc = 1:3
                    params.encPeriod = encPeriods{enc};
                    
                    % make save directory
                    pathExt = sprintf('t-%d_l-%d_cv-%d_enc-%d',t,l,cv,enc);
                    dataDir = fullfile(dataPath,pathExt);
                    saveDir = fullfile(basePath,pathExt);
                    if ~exist(saveDir,'dir')
                        mkdir(saveDir);
                    end
                    
                    % lastly, loop over phase, power, or powre and phase
                    % together (0=power,1=phase,2=both)
                    usePhase = [0,1,2];
                    for p = 1:3
                        
                        % get fixed parameters                                                
                        subjFiles = fullfile(dataDir,['*_lasso_pow.mat']);
                        if p==1
                            subjFiles = fullfile(dataDir,['*_lasso_phase.mat']);
                        elseif p==2
                            subjFiles = fullfile(dataDir,['*_lasso_powphase.mat']);
                        end           
                        subjFiles = dir(subjFiles);
                        
                        tBest   = NaN(1,length(subjFiles));
                        c       = NaN(1,length(subjFiles));
                        encBest = NaN(1,length(subjFiles));
                        aucs    = NaN(1,length(subjFiles));
                        
                        for f = 1:length(subjFiles)
                          subjData = load(fullfile(dataDir,subjFiles(f).name));
                          c(f) = mode([subjData.res.lambda{:}]);
                          tBest(f) = mode([subjData.res.tBest{:}]);
                          encBest(f) = mode([subjData.res.encBest{:}]);
                          aucs(f) = subjData.AUC;
                        end
                        
                        
                        include = cellfun(@isempty,strfind({subjFiles.name},subjs{s}));
                        params.timeBins  = timeBins{t}(mode(tBest(include)),:);
                        params.encPeriod = encPeriods{mode(encBest(include))};
                        params.Cs        = mode(c(include));
                        
                        params.usePhase = usePhase(p);
                        YC1_watrous_phase(subjs{s},params,saveDir);
                    end
                end
            end
        end
    end
end





