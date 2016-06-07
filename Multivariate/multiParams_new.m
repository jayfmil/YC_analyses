function params = multiParams_new()

% frequency bins to use
params.freqBins = [1 3;3 9;40 70;70 200];

timeBins = {[1 1000;1001 4000;4001 5000]};
fun = @(x) ceil((x+1000)/20);
timeBins = cellfun(fun,timeBins,'uniformoutput',false);
params.timeBins = timeBins{1};

% regions. If empty, use all electrodes. Choices right now are:
%          'mtl', 'hipp', 'ec', 'all', 'ca1','ca3','dg','sub','phc','prc',
%          'frontal','parietal','occipital','temporal','limbic'
params.region = '';


% filter to events of interest
params.eventFilter = @(events)allEncodingEvents(events);
params.basePath    = '/scratch/jfm2/YC1/multi';

% save out binned power to mat file?
params.savePower = 1;

% path to power data
params.powerPath = '/data10/scratch/jfm2/power50freqs';

% In YC1/2, each item has two encoding periods. Which to choose from?
params.encPeriods = [1,2];

% number of cross validation folds
params.useKfold = false;
params.k = 10;

% cross validation field ('blocknum' is really individiual item, 'block' is
% two consecutive items).
params.cvField = 'blocknum';

% nested cross validation field
params.nestedCvField = 'blocknum';

% what type of normalization
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

params.diffRegress = 1;

params.Gs = [];
params.usePhase = 0;
params.useWatrous = 0;

% permute the Y, usually in the process of creating a chance distribution
params.doPermute = 0;

function eventMask = allEncodingEvents(events)
eventMask = (strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD'));

function eventMask = UpperLowerThird_EncodingEvents(events)

% filter to just encoding events
allEncoding = (strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD'));

% determine thresholds for upper third of performance and lower third
perfErr = [events(allEncoding).testError];
lowThresh = prctile(perfErr,33);
highThresh = prctile(perfErr,67);

% restrict to just upper or lower third events
eventMask = (perfErr <= lowThresh) | (perfErr >= highThresh);
