function params = multiParams()

% Use native matlab functions or use liblinear (hint: liblinear is much
% faster). 'liblinear' or 'matlab'
params.library = 'liblinear';

% frequency bins to use
params.freqBins = [1 3;3 9;40 70;70 200];

timeBins = {[1 1000;1001 2000;2001 3000;3001 4000;4001 5000]};
fun = @(x) ceil((x+1000)/20);
params.timeBins = cellfun(fun,timeBins,'uniformoutput',false);

% for the reporting functions (makeSubjectReports and weightsByRegion),
% timeBinLabels is used to make label specific timepoints. If you don't
% want the labels, make it ''. If not empty, must be the same size as
% params.timeBins.
params.timeBinLabels = {'Pre','Spin','Drive1','Drive2','Drive3','Wait'};

% regions. If empty, use all electrodes. Choices right now are:
%          'mtl', 'hipp', 'ec', 'all', 'ca1','ca3','dg','sub','phc','prc',
%          'frontal','parietal','occipital','temporal','limbic'
params.region = '';

% individual model for each time bin, or all in one model
params.modelEachTime = 1;

% filter to events of interest
params.eventFilter = @(events)allEncodingEvents(events);
params.basePath    = '/data10/scratch/jfm2/YC1/multi';

% save out binned power to mat file?
params.savePower = 1;

% do binary classification or continuous regression
params.doBinary = 1;

% use original power (0) data or use regression corrected (1). 
params.useCorrectedPower = 0;

% path to power data
params.powerPath = '/data10/scratch/jfm2/power50freqs';

% In YC1/2, each item has two encoding periods. When selecting features, do
% we use just the first period ('first'), just the second period
% ('second'), use both but keep them as seperate observations ('both'), or
% use both but combine them ('combine'), or use the average ('average')
params.encPeriod = 'second';

% cross validation strictness.
%   0 = use all the data to calculate optimal lambda, and apply the model
%       with this lamda to all cross val folds
%   1 = Calculate a new lamda for each cross val. This is technically
%       more correct because then the test data is totall removed from the
%       train data (i.e., no peaking when calculating optimal penalty).
%       This takes a lot longer.
params.crossValStrictness = 1;

% number of cross validation folds
params.useKfold = false;
params.percentCV = 4/48;

% cross validation field ('blocknum' is really individiual item, 'block' is
% two consecutive items).
params.cvField = 'blocknum';

% nested cross validation field
params.nestedCvField = 'blocknum';

params.auc_prctileThresh = 100;

% ADD MULTIPLE ROUNDS OF SUB SAMPLING
% params.subSampNum = 10

% use lasso (L1) normaliztion or L2 normaliztion? L1 zeros out most of the
% weights, L2 minimizes but doesn't zero out. If L1, uses lassoglm(). If
% L2, uses logRegFun().
params.normType = 'L1';

% if empty, lambda will be computed using the specified cross validation
% strictness and nCV
params.lambda = [];

% what lassoglm Alpha value to use? 1 = lasso (L1), 0 similar to ridge
% (L2). In between is elastic net. This is only applicable if normType ==
% 'L1'
params.alpha = 1;

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

params.Cs = [];
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
