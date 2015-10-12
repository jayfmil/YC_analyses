function YC1_powerByObjLocation(subjs,params)
%function YC1_subject_summary(subjs)
%
% Univariate analyses correlating power and performance in YC1.
%
% Inputs:   subj - cell array of subjects to use (if empty, use all)
%         params - parameters to use, most importantly, what frequencies
%                  and time bin to use (if empty, use univarParams()).
%
% Outputs: none, saves results to params.basePath


% use default params if none given
if ~exist('params','var') || isempty(params)
    params = univarParams();
end

if size(params.timeBins,1) > 1
    fprintf('Multiple time bins not supported\n')
    return
end

% save directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
region = params.region;
if isempty(region);region = 'all';end
saveDir = fullfile(params.basePath,f(params.useCorrectedPower,y),region);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% see if this was submitted with an open pool. If so, parallel on the level
% of subjects. Otherwise, will loop over subjects one by one.
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    parfor s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        YC1_univarStats(subjs{s},params,saveDir);
    end
else
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        YC1_univarStats(subjs{s},params,saveDir);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function YC1_univarStats(subj,params,saveDir)
% function YC1_univarStats(subj,params,saveDir)
%
% For a given subject, loads power for electrodes in the specified region.
% Bins into specified time and frequencies. Correlates power with
% performance in the task.

% load tal structure
try
    tal = getBipolarSubjElecs(subj,params.doBipol,1,params.excludeEpiElecs);
    tal = filterTalByRegion(tal,params.region);
    if isempty(tal)
        fprintf('No %s electrode for %s.\n',params.region,subj)
        return
    end
    
    % load power parameters
    powParams = load(fullfile(params.powerPath,'params.mat'));
    
    % Setting time bins for convenience:
    tEnds     = (powParams.params.pow.timeWin:powParams.params.pow.timeStep:powParams.params.eeg.durationMS)+powParams.params.eeg.offsetMS;
    tStarts   = tEnds - powParams.params.pow.timeWin+1;
    powParams.timeBins = [tStarts' tEnds'];
    
    % load events
    events = get_sub_events('RAM_YC1',subj);
    
    % add the test error to the learning trials
    events  = addErrorField(events);
    
    % filter to events of interest
    eventsToUse = params.eventFilter(events);
    
    % response error for each trial
    er = [events(eventsToUse).testError];
    
    % also bin response error by region of evironment
    objLocs = vertcat(events(eventsToUse).objLocs);
    edgesX=linspace(0,28.4,params.numXBins+1);
    edgesY=linspace(0,14,params.numYBins+1);
    [~,xBin] = histc(abs(objLocs(:,1)),edgesX);
    [~,yBin] = histc(abs(objLocs(:,2)),edgesY);
    
    
    thresh = median(er);
    cond1 = er < thresh;
    cond2 = er >= thresh;
    
    % see if we have enough events
    if sum(cond1) < 5
        fprintf('Only %d events for %s in cond1 using %s. Skipping subject.\n', sum(cond1),subj,func2str(ana_func))
        return
    end
    
    if sum(cond2) < 5
        fprintf('Only %d events for %s in cond2 %s. Skipping subject.\n', sum(cond2),subj,func2str(ana_func))
        return
    end
    
    % load the power for this region
    region = params.region;
    if isempty(region)
        region = 'all';
    end
    fprintf('Calculating average power for %d %s elecs.\n',length(tal),region)
    powerData = loadAllPower(tal,subj,events,params.freqBins,params.timeBins,powParams,eventsToUse,params);
    powerData = permute(powerData,[3 4 1 2]);
    nElecs = size(powerData,2);
    
    % initial structure to hold results
    res = cell2struct(cell(1,length(params.freqBinLabels)),params.freqBinLabels,2);
    for field = params.freqBinLabels
        res.(field{1}) = struct('r',NaN(1,nElecs),'pCorr',NaN(1,nElecs),...
            'tstat',NaN(1,nElecs),'sd',NaN(1,nElecs),...
            'df',NaN(1,nElecs),'p_ttest',NaN(1,nElecs),...
            'meanCond1',NaN(1,nElecs),'meanCond2',NaN(1,nElecs),...
            'diffMap',NaN(params.numXBins,params.numYBins,nElecs),...
            'diffMapCorr',NaN(params.numXBins,params.numYBins,nElecs));
    end
    % tagNames = cell(1,size(elecs,1));
    
    % loop over each electrode in region
    for e = 1:size(powerData,2)
        
        % and each frequency
        for f = 1:size(powerData,3)
            
            % field name to save in res structure
            field = params.freqBinLabels{f};
            
            % power for this elec and freq
            pow = powerData(:,e,f)';
            
            % get the mean error for each bin
%             ind = sub2ind([params.numEnvBins params.numEnvBins],xBin,yBin);
%             powByLoc = NaN(params.numEnvBins^2,2);
%             m = grpstats(pow(cond1),ind(cond1));
%             powByLoc(unique(ind(cond1)),1) = m;
%             m2 = grpstats(pow(cond2),ind(cond2));
%             powByLoc(unique(ind(cond2)),2) = m2;
%             diffMap = NaN(params.numEnvBins,params.numEnvBins);
%             diffMap(1:(params.numEnvBins^2)) = powByLoc(:,1) - powByLoc(:,2);
%             res.(field).diffMap(:,:,e) = diffMap;
            
            % ttest and correlation for objects by bin of environment
            ind = sub2ind([params.numXBins params.numYBins],xBin,yBin);
            diffMap = NaN(params.numXBins,params.numYBins);
            diffMapCorr = NaN(params.numXBins,params.numYBins);
            uniqInds = unique(ind);
            for bin = uniqInds'
                
                % corr
                [r] = corr(er(ind==bin)',pow(ind==bin)');
                diffMapCorr(ind(ind==bin)) = r;
                
                % ttest
                [~,~,~,s] = ttest2(pow(ind==bin & cond1'),pow(ind==bin & cond2'));
                diffMap(ind(ind==bin)) = s.tstat;
            end
            res.(field).diffMapCorr(:,:,e) = diffMapCorr;
            res.(field).diffMap(:,:,e) = diffMap;
            
            % correlation between power and performance
            bad = isnan(er) | isnan(pow);
            [res.(field).r(e),res.(field).pCorr(e)] = corr(er(~bad)', pow(~bad)');
            
            % ttest between power in each condition
            [~,p,~,s] = ttest2(pow(cond1),pow(cond2));
            res.(field).tstat(e) = s.tstat;
            res.(field).sd(e) = s.sd;
            res.(field).df(e) = s.df;
            res.(field).p_ttest(e) = p;
            
            % mean power in each condition
            res.(field).meanCond1 = nanmean(pow(cond1));
            res.(field).meanCond2 = nanmean(pow(cond2));
            
        end
    end
    
    % save it to file
    fname = fullfile(saveDir,[subj '.mat']);
    save(fname,'res','tal')
catch
    fprintf('Error processing %s.\n',subj)
    return
end


function powerData = loadAllPower(tal,subj,events,freqBins,timeBins,powParams,eventsToUse,params)

nFreqs = size(freqBins,1);
nTimes = size(timeBins,1);
nEvents = sum(eventsToUse);
nElecs = length(tal);
powerData = NaN(nFreqs,nTimes,nEvents,nElecs);

% when loading power, use either original power or power with effect of
% trial number removed.
powField = 'pow';
if params.useCorrectedPower
    powField = 'powCorr';
end
for e = 1:nElecs
    elecNum = tal(e).channel;
    
    basePath  = '/data10/scratch/jfm2/RAM/biomarker/power/';
    subjPath  = fullfile(basePath,subj);
    sessions = unique([events.session]);
    subjPow  = [];
    for s = 1:length(sessions)
        fname = fullfile(subjPath,'RAM_YC1_events',num2str(sessions(s)),[num2str(elecNum(1)),'-',num2str(elecNum(2)),'.mat']);
        sessPow = load(fname);
        subjPow = cat(3,subjPow,sessPow.sessOutput.(powField));
    end
    
    if length(eventsToUse) ~= size(subjPow,3)
        fprintf('Number of events does not match size of power matrix for %s!.\n',subj)
        return
    end
    subjPow = subjPow(:,:,eventsToUse);
    
    % average frequencies
    if nFreqs ~= length(powParams.params.pow.freqs)
        tmpPower = NaN(nFreqs,size(subjPow,2),size(subjPow,3));
        for f = 1:nFreqs
            fInds = powParams.params.pow.freqs >= freqBins(f,1) & powParams.params.pow.freqs < freqBins(f,2);
            tmpPower(f,:,:) = nanmean(subjPow(fInds,:,:),1);
        end
        subjPow = tmpPower;
    end
    
    % average times
    tmpPower = NaN(nFreqs,nTimes,size(subjPow,3));
    for t = 1:nTimes
        tInds = powParams.timeBins(:,1) >= timeBins(t,1) & powParams.timeBins(:,2) <= timeBins(t,2);
        tmpPower(:,t,:) = nanmean(subjPow(:,tInds,:),2);
    end
    powerData(:,:,:,e) = tmpPower;
end

function events = addErrorField(events)
% add testError field
% add inner field (1 = inner region, 0 = outer region)

testInd = strcmp({events.type},'NAV_TEST');
recEvents = events(testInd);
[events.testError] = deal(NaN);
[events.recalled] = deal(NaN);
[events.inner] = deal(NaN);
[events.innerOuterBand] = deal(NaN);
sessVec = [events.session];
trialVec = [events.blocknum];
for rec = 1:length(recEvents);
    session = recEvents(rec).session;
    trial = recEvents(rec).blocknum;
    err = recEvents(rec).respPerformanceFactor;
    ind = sessVec == session & trialVec == trial;
    [events(ind).testError] = deal(err);
    [events(ind).inner] = deal(abs(recEvents(rec).objLocs(1)) < 568/30 && abs(recEvents(rec).objLocs(2)) < 7);
    [events(ind).innerOuterBand] = deal(abs(recEvents(rec).objLocs(1)) < 14);
end





