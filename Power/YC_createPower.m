function YC_createPower()
% function YC_createPower(saveDir)
% credit: Youssef Ezzyat, minor modifications by Jonathan Miller
%
% Creates power for each electrode for each subject in the RAM_YC1
% experiment.
%
% Note and acknowledgements: This borrows heavily from Youssef Ezzyat's
% code located in svn at rhino.psych.upenn.edu/home/svn/matlab/ye_code. His
% code is a dependancy for this to run.

% Youssef's comments below:
% creates a parameters struct for use in YE's eeg code. this params struct
% gets passed to ComputePow, for input to gete_ms and multiphasevec3. the
% function returns a struct variable params that contains fields specifying
% parameters for raw eeg voltage extraction, power computation and binning.
% the function also saves this params variable in a params.mat file in the
% directory specified in the input 'SaveDir'. this directory should be the
% one where the data is saved that was created using the specified
% parameters.
%
%
%   params.eeg.durationMS   length of the window within which to compute power
%   params.eeg.offsetMS     length of time by which durationMS precedes/follows event onset
%   params.eeg.bufferMS     size of buffer placed around epoch of interest
%   params.eeg.filtfreq     specific frequency to filter (e.g. [58 62] for 60Hz noise)
%   params.eeg.filttype     filter type (default in gete_ms = 'stop')
%   params.eeg.filtorder    filter order (default in gete_ms = 1)
%   params.eeg.sampFreq     resampling frequency in Hz (e.g. 1000 or 500)
%
%   params.pow.freqs        frequencies at which to compute power (e.g. [2 4 8...] or logspace(log10(2), log10(200),50)
%   params.pow.logTrans     1/0 flag indicating whether to log transform power vals
%   params.pow.wavenum      number of wavelets (typical = 7)
%   params.pow.timeWin      window size over which to average power when
%                           binning data
%   params.pow.timeStep     step size to jump when binning data
%   params.pow.freqBins     bins within which to average when binning


%%%%% create YC1 power using wavelets
task = 'RAM_YC1';
subjs = get_subs(task);

params.doBipol = 1;
params.eeg.durationMS   = 8000;
params.eeg.offsetMS     = -1000;
params.eeg.bufferMS     = 2000;
params.eeg.filtfreq     = [58 62];
params.eeg.filttype     = 'stop';
params.eeg.filtorder    = 4;
params.eeg.sampFreq     = 500;
params.eeg.kurtThr      = 4;
% params.pow.freqs        = logspace(log10(3),log10(180),8);
params.pow.logTrans     = 1;
% params.pow.type         = 'wavelet';
params.pow.wavenum      = 5;
params.pow.timeWin      = 20;
params.pow.timeStep     = 20;
% params.pow.freqBins     = logspace(log10(3),log10(180),8);

params.useGetPhasePow   = 1;

freqSettings = 'fifty'
if strcmp(freqSettings,'eight')
    params.pow.freqs        = logspace(log10(3),log10(180),8);
    params.pow.freqBins     = logspace(log10(3),log10(180),8);
    params.pow.type         = 'wavelet';
    params.savedir          = '/data10/scratch/jfm2/power8freqs';
elseif strcmp(freqSettings,'eightNew')
    params.pow.freqs        = logspace(log10(1),log10(200),8);
    params.pow.freqBins     = logspace(log10(1),log10(200),8);
    params.pow.type         = 'wavelet';
    params.savedir          = '/data10/scratch/jfm2/power8freqs_new';
elseif strcmp(freqSettings,'sixteen')
    params.pow.freqs        = logspace(log10(1),log10(200),16);
    params.pow.freqBins     = logspace(log10(1),log10(200),16);
    params.pow.type         = 'wavelet';
    params.savedir          = '/data10/scratch/jfm2/power16freqs';    
elseif strcmp(freqSettings,'fifty')
    params.pow.freqs        = logspace(log10(1),log10(200),50);
    params.pow.freqBins     = logspace(log10(1),log10(200),50);
    params.pow.type         = 'wavelet';
    params.savedir          = '/data10/scratch/jfm2/power50freqs';     
elseif strcmp(freqSettings,'four')
    params.pow.freqs        = [1 3;3 12;40 70;70 200];
    params.pow.freqBins     = [1 3;3 12;40 70;70 200];
    params.pow.type         = 'hilbert';
    params.savedir          = '/data10/scratch/jfm2/power4bins_hilbert';     
end

if ~params.doBipol
    params.savedir = [params.savedir,'_avg_reref'];
end

% for fft slep power
params.pow.bandwidth    = 2;
params.pow.winSize      = .4;
params.pow.winStep      = .05;

% also save out power with the effect of trial number regressed away?
params.regressTrialNumber = 1;

% if YC1, use all learning trials to zscore
params.eventsYC1        = @(events)strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD');

% If YC2, use all non-stim learning trials to zscore
params.eventsYC2        = @(events)[events.isStim]==0&(strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD'));
% params.savedir          = '/data10/scratch/jfm2/power8freqs';
cd_mkdir(params.savedir); save(['params_',task,'.mat'],'params');

% compute powers
fileExt = '';
computePower(task,subjs,params,fileExt)

%%%%% create YC2 power using wavelets
% task = 'RAM_YC2';
% subjs = get_subs(task);

% compute powers
% fileExt = '';
% computePower(task,subjs,params,fileExt)
% cd_mkdir(params.savedir); save(['params_',task,'.mat'],'params');
% 
% %%%% create YC2 post stim power with fft_slep method
% params.pow.type         = 'fft_slep';
% params.eeg.offsetMS     = 5100;
% params.eeg.durationMS   = 1900;
% params.eeg.wavenum = NaN;
% params.eeg.timeWin = NaN;
% params.eeg.timeStep = NaN;
% 
% cd_mkdir(params.savedir); save(['params_post_',task,'.mat'],'params');
% 
% % compute powers
% fileExt = '_post';
% computePower(task,subjs,params,fileExt)
% 
% %%%% create YC2 pre stim power with fft_slep method
% params.pow.type         = 'fft_slep';
% params.eeg.offsetMS     = -2000;
% params.eeg.durationMS   = 1900;
% 
% 
% cd_mkdir(params.savedir); save(['params_pre_',task,'.mat'],'params');
% 
% % compute powers
% fileExt = '_pre';
% computePower(task,subjs,params,fileExt)