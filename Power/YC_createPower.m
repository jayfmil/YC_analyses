function YC_createPower(task,subjs)
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

if ~ismember(task,{'RAM_YC1','RAM_YC2'});
    fprintf('Supported tasks names are RAM_YC1 and RAM_YC2.\n')
    return
end

% if subjects not given, use all
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs(task);
end

params.eeg.durationMS   = 8000;
if strcmp(task,'RAM_YC2')
    params.eeg.durationMS   = 8000;
end
params.eeg.offsetMS     = -1000;
params.eeg.bufferMS     = 2000;
params.eeg.filtfreq     = [58 62];
params.eeg.filttype     = 'stop';
params.eeg.filtorder    = 1;
params.eeg.sampFreq     = 500;
params.eeg.kurtThr      = 4;
params.pow.freqs        = logspace(log10(1),log10(200),50);
params.pow.logTrans     = 1;
params.pow.type         = 'wavelet';
params.pow.wavenum      = 7;
params.pow.timeWin      = 100;
params.pow.timeStep     = 100;
params.pow.freqBins     = logspace(log10(1),log10(200),50);

% also save out power with the effect of trial number regressed away?
params.regressTrialNumber = 1;

% if YC1, use all learning trials to zscore
params.eventsYC1        = @(events)strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD');

% If YC2, use all non-stim learning trials to zscore
params.eventsYC2        = @(events)[events.isStim]==0&(strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD'));
params.savedir          = '/data10/scratch/jfm2/power';
cd_mkdir(params.savedir); save(['params_',task,'.mat'],'params');

% compute powers
computePower(task,subjs,params)






