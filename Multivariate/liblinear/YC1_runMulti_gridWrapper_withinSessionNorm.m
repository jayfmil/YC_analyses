function YC1_runMulti_gridWrapper_withinSessionNorm(subjs)

basePath = '/scratch/jfm2/YC1/multi/watrousGrid';
basePath = '/scratch/jfm2/YC1/multi/pow8_grid_Cs_iterate_withinSessionNorm'
powerPath = '/scratch/jfm2/power8freqs/';

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% use default params if none given
% if ~exist('params','var') || isempty(params)
%     params = multiParams();
% end
% params.powerPath = powerPath;
% params.crossValStrictness = 1;
% params.useKfold = 0;
% params.cvField  = 'session';
% params.useWatrous = 1;
% params.basePath = basePath;
% params.region = 'mtl';
% params.overwrite = 1;

Cs = {};
Cs{1} = logspace(log10(1e-2),log10(1e4),22);
Cs{2} = logspace(log10(1e-6),log10(1e4),22);

% loop over subjets
parfor s = 1:length(subjs)
    
    % loop over time bin settings
    %     timeBins = {[1 20;21 40;41 60;61 80;81 100],[1 20],[81 100],[1 100]};
        timeBins = {[1 1000;1001 2000;2001 3000;3001 4000;4001 5000],...
                    [1 1000],...
                    [4001 5000],...
                    [1 5000]};
        fun = @(x) ceil((x+1000)/20);
        timeBins = cellfun(fun,timeBins,'uniformoutput',false);
%     timeBins = {[[1:100]' [1:100]']};
    for t = 1:length(timeBins)
%         params.timeBins = timeBins{t};
        timeBin  = timeBins{t};
        
        % loop over regularization method
        normTypes = {'L1','L2'};
        for l = 1:2
%             params.normType = normTypes{l};
            normType = normTypes{l};
            
            % loop over inner cross val scheme
            nestedCvFields = {'session','blocknum'};
%             nestedCvFields = {'blocknum'};
            for cv = 1:2
%                 params.nestedCvField = nestedCvFields{cv};
                nestedCvField = nestedCvFields{cv};
                
                % loop over encoding period selection
                encPeriods = {'first','second','combined'};                
                for enc = 3:3
%                     params.encPeriod = encPeriods{enc};
                    encPeriod = encPeriods{enc};
                    
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
                        for lambdas = 1:length(Cs{l})
                            
                            params = multiParams();
                            params.powerPath = powerPath;
                            params.crossValStrictness = 1;
                            params.useKfold = 0;
                            params.cvField  = 'session';
%                             params.cvField  = 'blocknum';
                            params.useWatrous = 0;
                            params.basePath = basePath;
                            params.region = 'mtl';
                            params.region = '';
                            params.overwrite = 1;
                            params.timeBins = timeBin;
                            params.normType = normType;
                            params.nestedCvField = nestedCvField;
                            params.encPeriod = encPeriod;
                            
                            pathExt = sprintf('Cs_%d',lambdas);
                            saveDirSub = fullfile(saveDir,pathExt);
                            params.Cs = Cs{l}(lambdas);
                            
                            if ~exist(saveDirSub,'dir')
                                mkdir(saveDirSub);
                            end
                            params.usePhase = usePhase(p);
                            YC1_watrous_phase_withinSessionNorm(subjs{s},params,saveDirSub);
                        end
                    end
                end
            end
        end
    end
end







