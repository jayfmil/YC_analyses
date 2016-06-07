function params = TH_multiParams()



params.freqBins = [1 3;3 12;40 70;70 200];
params.timeBins = [51 125];

% regions. If empty, use all electrodes. Choices right now are:
%          'mtl', 'hipp', 'ec', 'all', 'ca1','ca3','dg','sub','phc','prc',
%          'frontal','parietal','occipital','temporal','limbic'
params.region = '';

% filter to events of interest
params.eventFilter = @(events)allEncodingEvents(events);
params.basePath    = '/scratch/jfm2/TH1/multi';

% save out binned power to mat file?
params.savePower = 1;

% do binary classification or continuous regression
params.doBinary = 1;

% use original power (0) data or use regression corrected (1). 
params.useCorrectedPower = 0;

% path to power data
params.powerPath = '/data10/scratch/jfm2/power50freqs';

% number of cross validation folds
params.useKfold = false;
params.k = 10;

% cross validation field ('blocknum' is really individiual item, 'block' is
% two consecutive items).
params.cvField = 'trial';

% nested cross validation field
params.nestedCvField = 'trial';

params.normType = 'L2';


% what regularization parameters to use
params.Cs = [];
if isempty(params.Cs)    
    if strcmpi(params.normType,'L1')
        params.Cs = logspace(log10(1e-2),log10(1e4),22);
    elseif strcmpi(params.normType,'L2')
        params.Cs = logspace(log10(1e-6),log10(1e4),22);
    end
end

% save the output to a file? Might not want to in some casesm for example
% creating a chance distribution
params.saveOutput = 1;

% default is to do nothing if the file already exists and saveOutout is
% true
params.overwrite = 0;

% load prebinned power
params.loadPower = 0;

% exclude epileptic electrodes?
params.excludeEpiElecs = 0;

params.correctThresh = NaN;

params.Gs = [];
params.usePhase = 0;
params.useWatrous = 0;

% permute the Y, usually in the process of creating a chance distribution
params.doPermute = 0;

function eventMask = allEncodingEvents(events)
eventMask = strcmp({events.type},'CHEST')  & ~cellfun('isempty',{events.item});
