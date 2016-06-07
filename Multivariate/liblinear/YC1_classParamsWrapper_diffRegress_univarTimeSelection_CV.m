function YC1_classParamsWrapper_diffRegress_univarTimeSelection_CV(useWatrous)

% get initial params
baseParams = multiParams;
baseParams.basePath = '/data10/scratch/jfm2/YC1/multi/acrossTrial_modelCompare_diffReg_univarTime_CV';
regions    = {'all','mtl'};
regions    = {'all'}
subjs = get_subs('RAM_YC1');
subjs = subjs(~strcmp(subjs,'R1065J'));
subjs = subjs(~strcmp(subjs,'R1047D'));
subjs = subjs(~strcmp(subjs,'R1062J'));
subjs = subjs(~strcmp(subjs,'R1096E'));
subjs = subjs(~strcmp(subjs,'R1037D'));
subjs = subjs(~strcmp(subjs,'R1051J'));

if ~useWatrous
    freqBins = {[1 12;40 200],[1 3;3 12;40 70;70 200],[]};
    powerPaths = {'/scratch/jfm2/power50freqs','/scratch/jfm2/power50freqs','/scratch/jfm2/power8freqs/'};
    timeBin = {[1 1000;1001 2000;2001 3000;3001 4000;4001 5000;5001 6000]};
    fun = @(x) ceil((x+1000)/20);
    timeBin = cellfun(fun,timeBin,'uniformoutput',false);
    regTypes = {'L1','L2','svm'};
    regTypes = {'L1','L2'};
    watrousData = 'nonWatrousData';
else
    timeBin  = [1 100];
    freqBins = {[]};
    regTypes = {'L1','L2','svm'};
    watrousData = 'WatrousData';
end

nFreqBins = length(freqBins);
nRegTypes = length(regTypes);

aucs = {};
classType = {};

for f = 1:nFreqBins
    params  = baseParams;
    freqBin = freqBins{f};
    params.powerPath = powerPaths{f};
    params.freqBins = freqBin;
    params.timeBins = timeBin{1};
    params.useWatrous = useWatrous;
    params.cvField = 'blocknum';
    params.diffRegress = 1;
    params.encPeriod = 'both'
    
%     params.powerPath = '/scratch/jfm2/power50freqs';
    
    for r = 1:length(regions)
        region = regions{r};
        params.region = region;
        
        for n = 1:nRegTypes
            regType = regTypes{n};
            params.normType = regType;
            
            pathExt = sprintf('%s_%dFreqs_region_%s_%s',regType,size(freqBin,1),region,watrousData);
            classType{end+1} = sprintf('%s_%dFreqs_region_%s',regType,size(freqBin,1),region);
            params.basePath = fullfile(baseParams.basePath,pathExt);            
            
            YC1_runMulti_phase_univarTimeSelection(subjs,params)
            aucs{end+1} = YC1_plotClassRes([],params,1);
            
        end
    end
end
keyboard