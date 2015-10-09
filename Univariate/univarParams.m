function params = univarParams()

% frequency bins to use
params.freqBins = [1 3;3 9;40 70;70 200];

% time points
params.timeBins = [1 5000];

% do bipolar
params.doBipol = 1;

% regions. If empty, use all electrodes. Choices right now are:
%          'mtl', 'hipp', 'ec', 'all', 'ca1','ca3','dg','sub','phc','prc',
%          'frontal','parietal','occipital','temporal','limbic'
params.region = '';

% filter to events of interest
params.eventFilter = @(events)allEncodingEvents(events);
params.basePath    = '/data10/scratch/jfm2/YC1/uni';

% save out binned power to mat file?
params.savePower = 1;

% use original power (0) data or use regression corrected (1). 
params.useCorrectedPower = 0;

% path to power data
params.powerPath = '/data10/scratch/jfm2/power';

% save the output to a file? Might not want to in some casesm for example
% creating a chance distribution
params.saveOutput = 1;

% default is to do nothing if the file already exists and saveOutout is
% true
params.overwrite = 0;

% load prebinned power
params.loadPower = 0;

% exclude epileptic electrodes?
params.excludeEpiElecs = 1;

% permute the Y, usually in the process of creating a chance distribution
params.doPermute = 0;

function eventMask = allEncodingEvents(events)
eventMask = (strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD'));

function eventsMask = correctFilter(events, isCorrect)     

encInd = strcmp({events.type},'NAV_LEARN') | strcmp({events.type}, ...
                                                  'NAV_LEARN_HARD');
if ~isCorrect
  eventsMask = [events.recalled]==0 & encInd;
else
  eventsMask = [events.recalled]==1 & encInd;
end

function eventMask = UpperLowerThird_EncodingEvents(events)

% filter to just encoding events
allEncoding = (strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD'));

% determine thresholds for upper third of performance and lower third
perfErr = [events(allEncoding).testError];
lowThresh = prctile(perfErr,33);
highThresh = prctile(perfErr,67);

% restrict to just upper or lower third events
eventMask = (perfErr <= lowThresh) | (perfErr >= highThresh);
