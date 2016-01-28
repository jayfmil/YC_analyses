function PowMat = computeFFT_slepPow(EEG,params)

% some parameters:
params.pow.freqs        = logspace(log10(1),log10(200),50);
params.eeg.sampFreq     = 500;
params.pow.bandwidth    = 2;
params.pow.winSize      = .4;


winSize = params.pow.winSize*params.eeg.sampFreq; % winSize in samples

% figure out the length to make the fft, to get the freq resolution
% you need to return the freqs specified in params.pow.freqs
nfft = max(2^nextpow2(params.eeg.sampFreq/min(diff(params.pow.freqs))),winSize);
%nfft = max(2^nextpow2(winSize),winSize);


% don't include any buffer in the fft data
buffToRem = params.eeg.bufferMS*params.eeg.sampFreq/1000;
buffInd = [1:buffToRem, (size(EEG,2)-buffToRem+1):size(EEG,2)];
EEG(:,buffInd) = [];


% loop through the signal to compute FFT in windows of size params.pow.winSize
SignalLen = size(EEG,2);

% Compute the positions of the windows
if isnan(params.pow.winStep) % do every sample
    winStarts = 1:((SignalLen-winSize)+1);
    winEnds = winSize:SignalLen;
else % use a timestep
    winStep = params.pow.winStep*params.eeg.sampFreq;
    winStarts = 1:winStep:((SignalLen-winSize)+1);
    winEnds = winSize:winStep:SignalLen;
end

nEvents = size(EEG,1);

% define the sleppian windows
sWins = dpss(winSize, params.pow.bandwidth, (2*params.pow.bandwidth)-1); %sWins has dimensions (N x numSWins)
numSWins = size(sWins, 2);
sWins = sWins'; %now sWins has dimensions (numSWins x numSamples)

sWins    = sWins(:,:,ones(1,nEvents));
sWins    = permute(sWins,[3 2 1]);

% allocate space for fft, define freqs & nyquist
%nFreqs = length(params.pow.freqs);
Fs = 0:params.eeg.sampFreq/nfft:params.eeg.sampFreq;
Fs = Fs(1:nfft);
NyquistInds = Fs<=floor(params.eeg.sampFreq/2);
% get the fft freq index that is closest to our desired freq index,
% and also the one prior and after. take the average across all 3
% to get estimate of desired frequency
freqInds = nan(length(params.pow.freqs),1);
for iFreq = 1:length(params.pow.freqs)
    [v,c] = min(abs(Fs-params.pow.freqs(iFreq)));
    freqInds(iFreq) = c;
end

% initialize PowMat
PowMat = nan(length(params.pow.freqs), length(winStarts), nEvents);

% expand EEG to 3D for ease of multiplication w/slep wins
EEG = EEG(:,:,ones(1,numSWins));

for iInd = 1:length(winStarts)
    StartInd = winStarts(iInd);
    EndInd = winEnds(iInd);
    
    % multiply signal by windows
    tmpEEG = EEG(:,StartInd:EndInd,:);
    tmpEEGProj = tmpEEG.*sWins;
    
    % compute the fft
    FT = fft(tmpEEGProj,nfft,2)/nfft;
    FT = abs(FT).^2;
    FT = mean(FT,3);
    
    %FT = FT(:,freqInds)'; % select out freqs closest to those requested. also get prior and next freq, and avg the 3 together
    FT = cat(3,FT(:,freqInds),FT(:,freqInds-1),FT(:,freqInds+1));
    FT = nanmean(FT,3)';
    
    PowMat(:,iInd,:) = FT;
end

if params.pow.logTrans==1
    PowMat = log10(PowMat);
end



