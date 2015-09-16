function YC_plotElecTrialPower(subj,task,elecNums)

% parameters obviously
params.eeg.durationMS   = 2000;
params.eeg.offsetMS     = -1000;
params.eeg.bufferMS     = 2000;
params.eeg.filtfreq     = .1;
params.eeg.filtfreq     = [58 62];
params.eeg.filttype     = 'high';
params.eeg.filttype     = 'stop';
params.eeg.filtorder    = 1;
params.eeg.sampFreq     = 500;
params.eeg.kurtThr      = 4;
params.pow.freqs        = logspace(log10(1),log10(256),80);
params.pow.logTrans     = 1;
params.pow.type         = 'wavelet';
params.pow.wavenum      = 7;
params.pow.timeWin      = 100;
params.pow.timeStep     = 100;
params.pow.freqBins     = logspace(log10(1),log10(256),80);
params.events           = @(events)strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD');
% params.savedir          = saveDir;

% get events for subject
events   = get_sub_events(task,subj);
%events   = addErrorField(events);
sessions = unique([events.session]);
for sess = 1:length(sessions)
   

    % filter to session events
    sess_events = events([events.session]==sessions(sess));
    
    % filter to just learning trials
    sess_events = sess_events(strcmp({sess_events.type},'WORD'));
    
    % filter to non-stim
    sess_events = sess_events([sess_events.isStim]==0);
    
%    sess_events = sess_events([sess_events.serialpos]==1);

    % calculate power for this session and electrode
    eeg = ComputeEEG(elecNums,sess_events,params);
    pow = ComputePow(eeg,params);
    keyboard
    [~,pow1] = getphasepow(elecNums(1),sess_events,params.eeg.durationMS,params.eeg.offsetMS,params.eeg.bufferMS,'freqs',params.pow.freqs,'width',params.pow.wavenum);
    [~,pow2] = getphasepow(elecNums(2),sess_events,params.eeg.durationMS,params.eeg.offsetMS,params.eeg.bufferMS,'freqs',params.pow.freqs,'width',params.pow.wavenum);
    pow1 = permute(pow1,[2 3 1]);
    pow1(pow1<0) = eps;
    pow1 = log10(pow1);
    pow2 = permute(pow2,[2 3 1]);
    pow2(pow2<0) = eps;
    pow2 = log10(pow2);

    figure(1)
    clf
    cmap      = zeros(size(pow,3),3);
    cmap(:,1) = linspace(1,0,size(pow,3));
    for c = 1:size(cmap,1)
        hold on
        plot(nanmean(pow(:,:,c),2),'color',cmap(c,:));                
    end   
    f2 = 2.^(0:ceil(sqrt(max(params.pow.freqs))));
    f2 = f2(f2<=(max(round(params.pow.freqs))));
    inds = interp1(log10(logspace(log10(1),log10(max(f2)),length(params.pow.freqs))),1:length(params.pow.freqs),log10(f2),'linear');
    set(gca,'xtick',inds);
    set(gca,'xticklabel',f2);
    grid on
    set(gca,'gridlinestyle',':')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('log(power)','fontsize',16)
    set(gca,'fontsize',16)
    
%    recalled = [sess_events.testError] < median([sess_events.testError]);

    recalled = [sess_events.recalled] == 1;
    hfaPow   = nanmean(squeeze(nanmean(pow(params.pow.freqs>70&params.pow.freqs<200,:,:),2)),1);


    
    figure(2)
    clf
    plot(find(recalled),hfaPow(recalled),'.r','markersize',40)
    hold on
    plot(find(~recalled),hfaPow(~recalled),'.b','markersize',40)
    grid on
    set(gca,'gridlinestyle',':')
    xlabel('Learning Trial Number','fontsize',16)
    ylabel('HFA log(power)','fontsize',16)
    set(gca,'fontsize',16)
    keyboard

    figure(3)
    clf
    cmap      = zeros(size(pow1,3),3);
    cmap(:,1) = linspace(1,0,size(pow1,3));
    for c = 1:size(cmap,1)
        hold on
        plot(nanmean(pow1(:,:,c),2),'color',cmap(c,:));
    end
    f2 = 2.^(0:ceil(sqrt(max(params.pow.freqs))));
    f2 = f2(f2<=(max(round(params.pow.freqs))));
    inds = interp1(log10(logspace(log10(1),log10(max(f2)), ...
                                  length(params.pow.freqs))),1:length(params.pow.freqs),log10(f2),'linear');
    set(gca,'xtick',inds);
    set(gca,'xticklabel',f2);
    grid on
    set(gca,'gridlinestyle',':')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('log(power)','fontsize',16)
    set(gca,'fontsize',16)

    recalled = [sess_events.testError] < ...
        median([sess_events.testError]);
    hfaPow   = nanmean(squeeze(nanmean(pow1(params.pow.freqs>70& ...
                                           params.pow.freqs<200,:,:),2)),1);

    figure(4)
    clf
    plot(find(recalled),hfaPow(recalled),'.r','markersize',25)
    hold on
    plot(find(~recalled),hfaPow(~recalled),'.b','markersize',25)
    grid on
    set(gca,'gridlinestyle',':')
    xlabel('Learning Trial Number','fontsize',16)
    ylabel('HFA log(power)','fontsize',16)
    set(gca,'fontsize',16)
        

    figure(5)
    clf
    cmap      = zeros(size(pow2,3),3);
    cmap(:,1) = linspace(1,0,size(pow2,3));
    for c = 1:size(cmap,1)
        hold on
        plot(nanmean(pow2(:,:,c),2),'color',cmap(c,:));
    end
    f2 = 2.^(0:ceil(sqrt(max(params.pow.freqs))));
    f2 = f2(f2<=(max(round(params.pow.freqs))));
    inds = interp1(log10(logspace(log10(1),log10(max(f2)), ...
                                  length(params.pow.freqs))),1:length(params.pow.freqs),log10(f2),'linear');
    set(gca,'xtick',inds);
    set(gca,'xticklabel',f2);
    grid on
    set(gca,'gridlinestyle',':')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('log(power)','fontsize',16)
    set(gca,'fontsize',16)

    recalled = [sess_events.testError] < ...
        median([sess_events.testError]);
    hfaPow   = nanmean(squeeze(nanmean(pow2(params.pow.freqs>70& ...
                                           params.pow.freqs<200,:,:),2)),1);

    figure(6)
    clf
    plot(find(recalled),hfaPow(recalled),'.r','markersize',25)
    hold on
    plot(find(~recalled),hfaPow(~recalled),'.b','markersize',25)
    grid on
    set(gca,'gridlinestyle',':')
    xlabel('Learning Trial Number','fontsize',16)
    ylabel('HFA log(power)','fontsize',16)
    set(gca,'fontsize',16)

    keyboard
end

function events = addErrorField(events)
% add testError field
% add inner field (1 = inner region, 0 = outer region)

testInd = strcmp({events.type},'NAV_TEST');
recEvents = events(testInd);
[events.testError] = deal(NaN);
[events.recalled] = deal(NaN);
[events.inner] = deal(NaN);
sessVec = [events.session];
trialVec = [events.blocknum];
for rec = 1:length(recEvents);
    session = recEvents(rec).session;
    trial = recEvents(rec).blocknum;
    err = recEvents(rec).respPerformanceFactor;
    ind = sessVec == session & trialVec == trial;
    [events(ind).testError] = deal(err);
    [events(ind).inner] = deal(abs(recEvents(rec).objLocs(1)) < 568/30 && abs(recEvents(rec).objLocs(2)) < 7);
end









