function TH_classParamsWrapper(useWatrous)

% get initial params
baseParams = TH_multiParams;
regions    = {'all','mtl'};

if ~useWatrous
    freqBins = {[1 12;40 200],[1 3;3 12;40 70;70 200]};
    timeBin  = [51 125];
    regTypes = {'L1','L2','svm'};
    watrousData = 'nonWatrousData';
else
    timeBin  = [1 21];
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
    params.freqBins = freqBin;
    params.timeBins = timeBin;
    params.useWatrous = useWatrous;
    params.powerPath = '/scratch/jfm2/power50freqs';
    
    for r = 1:length(regions)
        region = regions{r};
        params.region = region;
        
        for n = 1:nRegTypes
            regType = regTypes{n};
            params.normType = regType;
            
            pathExt = sprintf('%s_%dFreqs_region_%s_%s',regType,size(freqBin,1),region,watrousData);
            classType{end+1} = sprintf('%s_%dFreqs_region_%s',regType,size(freqBin,1),region);
            params.basePath = fullfile(baseParams.basePath,pathExt);
            
            TH_runMulti([],params)
            aucs{end+1} = TH_plotClassRes([],params,1);
            
        end
    end
end
keyboard