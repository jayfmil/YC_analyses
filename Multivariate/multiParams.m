function params = multiParams()

% frequency bins to use
params.freqBins = [1 3;3 9;40 70;70 200];
% params.freqBins = [1 4;70 200];


% time bins to use (this depends on the RAM auto computed power time window)
params.timeBins = [-999 0;1 1000;1001 4000;4001 5000];

% filter to events of interest
params.eventFilter = @(events)allEncodingEvents(events);
% params.eventFilter = @(events)UpperLowerThird_EncodingEvents(events);

% save out binned power to mat file?
params.savePower = 1;

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
