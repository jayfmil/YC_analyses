function params = multiParams()

% frequency bins to use
params.freqBins = [1 3;3 9;40 70;70 200];
% params.freqBins = [1 9;40 200];

% time bins to use (this depends on the RAM auto computed power time window)
% params.timeBins = [-999 0;1 1000;1001 4000;4001 5000];
timeStep = 1000;
params.timeBins = [[-999:timeStep:5000]' [(-999+timeStep-1):timeStep:5000]'];
params.timeBins = [1 1000;1001 4000;4001 5000];

% individual model for each time bin, or all in one model
params.modelEachTime = 0;

% filter to events of interest
params.eventFilter = @(events)allEncodingEvents(events);
params.basePath    = '/data10/scratch/jfm2/YC1/multi/lassoReg_allEncoding_binary';

% save out binned power to mat file?
params.savePower = 1;

% do binary classification or continuous regression
params.doBinary = 1;

% use original power (0) data or use regression corrected (1)
params.useCorrectedPower = 1;
params.regressDir = '/data10/scratch/jfm2/YC1/multi/power/regress';
params.lassDir    = '/data10/scratch/jfm2/YC1/multi';

% cross validation strictness.
%   0 = use all the data to calculate optimal lambda, and apply the model
%       with this lamda to all cross val folds
%   1 = Calculate a new lamda for each cross val. This *might* technically
%       be more correct because then the test data is totall removed from
%       the train data. This takes a lot longer.
params.crossValStrictness = 0;

% save the output to a file? Might not want to in some casesm for example
% creating a chance distribution
params.saveOutput = 1;

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
