function YC_calcTrialPower(subjs)
%function YC_calcTrialPower(subjs)
%
% Calculate average power for each encording and test period for each
% electrode
%

saveDir = '/scratch/jfm2/YC1/averageTrialPower';
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end


params = [];
params.duration = 5000;
params.buffer = 2000;
params.offset = 0;
params.width = 7;
params.freqs = logspace(log10(1),log10(200),50);
params.doBipol = 1;

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
        YC1_trialPower(subjs{s},params,saveDir);
    end
else
    for s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        YC1_trialPower(subjs{s},params,saveDir);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function YC1_trialPower(subj,params,saveDir)


saveDir = fullfile(saveDir,subj);
if ~exist(saveDir,'dir')
  mkdir(saveDir)
end

% load tal structure
try
    
    % get all electrodes
    tal = getBipolarSubjElecs(subj,params.doBipol,1);
    
    
    % load events
    events    = get_sub_events('RAM_YC1',subj);
    testInds = strcmp({events.type},'NAV_TEST') | strcmp({events.type},'NAV_PRACTICE_TEST');
    

    % want to save out some event specific info
    % 'objLocs','respLocs','startLocs','respErrs'
    objLocs  = vertcat(events.objLocs);
    respLocs = vertcat(events.respLocs);
    startLocs = vertcat(events.startLocs);

    duration = NaN(1,length(events));
    duration(~testInds) = params.duration;
    for e = 1:length(tal)
        
        % holds events x freqs for this elec
        powAll = NaN(length(events),length(params.freqs));
        
        % get power for all learn events (same duration)
        chan = tal(e).channel;
        [~,pow] = getphasepow_bipol(chan,events(~testInds),params.duration,params.offset,params.buffer,'width',params.width','freqs',params.freqs);
        pow(pow<=0) = eps;
        pow = log10(pow);
        pow = nanmean(pow,3);
        powAll(~testInds,:) = pow;
        
        % now get power for each test event
        testEvs = events(testInds);
        inds = find(testInds);
        for testEv = 1:sum(testInds)
            thisDur = round((testEvs(testEv).respReactionTime + testEvs(testEv).respTravelTime) * 1000);
            duration(inds(testEv)) = thisDur;
            [~,pow] = getphasepow_bipol(chan,testEvs(testEv),thisDur,params.offset,params.buffer,'width',params.width','freqs',params.freqs);
            pow(pow<=0) = eps;
            pow = log10(pow);
            pow = nanmean(pow,3);
            powAll(inds(testEv),:) = pow;            
        end
        
        % zscore by session
        sessions = [events.session];
        zPowAll = zscore_local(powAll,sessions);
        
        
        % save it to file
        if length(chan)>1
          fname = [num2str(chan(1)) '-' num2str(chan(2)) '.mat'];
        else
          fname = [num2str(chan) '.mat'];
        end      
        fname = fullfile(saveDir,fname);
        save(fname,'powAll','zPowAll','objLocs','respLocs','startLocs','respErrs','sessions')
    end
catch e
  keyboard
  fprintf('Error processing %s.\n',subj)
  return
end


function zpow = zscore_local(pow,sessions)


zpow = NaN(size(pow));
uniqSessions = unique(sessions);
for s = 1:length(uniqSessions)
  
  sessInds = sessions==uniqSessions(s);
  powSess  = pow(sessInds,:);
  m = nanmean(powSess,1);
  s = nanstd(powSess);

  zPowSess = (powSess - repmat(m,[sum(sessInds),1])) ./ repmat(s,[sum(sessInds),1]);
  zpow(sessInds,:) = zPowSess;

end


