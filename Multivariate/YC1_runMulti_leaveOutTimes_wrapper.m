function YC1_runMulti_leaveOutTimes_wrapper
%
% Run classification including many time bins. Then leave out each time bin
% and re-run classification with the bin excluded. The goal is to see how
% the removal of specific time bins affects classification performance.
%

if isunix && ~ismac
    open_rhino2_pool(50,'12G');
end

% get basic parameters
params = multiParams();
params.modelEachTime = 0;
basePath = '/data10/scratch/jfm2/YC1/multi/excludeTimes';
timeBins = [-999     0;...
            1     1000;...
            1001  2000;...
            2001  3000;...
            3001  4000;...
            4001  5000;...
            5001  6000]; 
timeBinLabels = {'Pre', 'Spin', 'Drive1', 'Drive2', 'Drive3', 'Wait', 'Post'};

% get subjects and exclude R1001P because subject takes FOREVER to run
subjs = get_subs('RAM_YC1');
subjs = subjs(~strcmp(subjs,'R1001P'));

% First model with all time bins included
params.timeBins = timeBins;
params.timeBinLabels = timeBinLabels;
params.basePath = fullfile(basePath,'allTimes');
YC1_runMulti(subjs,params);
YC1_createChanceDists(subjs,params)

% loop over time bins, remove from classification
timesToRemove = {'Pre', 'Spin', 'Drive1', 'Drive2', 'Drive3', 'Wait', 'Post','Drive'};
for t = 1:length(timesToRemove)
    
    % Note: I'm using strfind to match all the timeBins that start with a
    % certain string so that I can catch all the Drive* at once. Careful
    % with overlapping strings
    binsToKeep = cellfun('isempty',strfind(timeBinLabels,timesToRemove{t}));
    params.timeBins = timeBins(binsToKeep,:);
    params.timeBinLabels = timeBinLabels(binsToKeep);
    pathStr = ['allTime_exc',timesToRemove{t}];
    params.basePath = fullfile(basePath,pathStr);
    YC1_runMulti(subjs,params);
    YC1_createChanceDists(subjs,params)
end
