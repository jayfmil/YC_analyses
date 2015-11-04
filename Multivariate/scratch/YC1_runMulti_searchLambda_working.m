function YC1_runMulti_searchLambda_working
%
%
%

if isunix && ~ismac
    p = gcp('nocreate');
    if isempty(p)
        open_rhino2_pool(50,'12G');
    end
end


% get subjects
subjs  = get_subs('RAM_YC1');

parfor s = 1:length(subjs)
    try
        % possible lambda value
        lambdas = logspace(log10(.01),log10(1),25);
        
        % get basic parameters
        params = multiParams();
        params.modelEachTime = 0;
        
        timeStep = 250;
        params.timeBins = [[-999:timeStep:0]' [(-999+timeStep-1):timeStep:0]'];
        params.timeBinLabels = strread(num2str(1:size(params.timeBins,1)),'%s')';
        params.freqBins = [[1:50]' [1:50]'];
        params.nCV = 10;
        params.basePath = '/data10/scratch/jfm2/YC1/multi/lambdaSearchPre50';
        params.saveOutput = 0;
        
        % save directory
        f = @(x,y) y{double(x)+1};
        y = {'OrigPower','CorrectedPower'};
        saveDir = fullfile(params.basePath,f(params.useCorrectedPower,y));
        if ~exist(saveDir,'dir')
            mkdir(saveDir);
        end
        
        aucs = NaN(1,length(lambdas));
        for l = 1:length(lambdas)
            
            if l == 1
                params.savePower = 1;
                params.loadPower = 0;
            else
                params.savePower = 0;
                params.loadPower = 1;
            end
            
            params.lambda = lambdas(l);
            [~,aucs(l)] = YC1_runMulti_subj(subjs{s},params,saveDir)
            
        end
        
        % find lambda with best performance, rerun and save output
        [~,ind]           = max(aucs);
        lambda            = lambdas(ind);
        params.lambda     = lambda;
        params.saveOutput = 1
        
        YC1_runMulti_subj(subjs{s},params,saveDir);
    end
end