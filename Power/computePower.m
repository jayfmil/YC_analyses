function computePower(task,subjs,params,fileExt)
% function computePower(task,subjs,params)
% credit: Youssef Ezzyat, minor modifications by Jonathan Miller. The
% largest change is the ability to regress out trial number from the power.

% loop over each subj
for s = 1:length(subjs)
    
    % stupid try because getBipolarSubjElecs fails if not file found
    try
           
        subj = subjs{s};
        events = get_sub_events(task,subj);        
        sessions = unique([events.session]);
        tal = getBipolarSubjElecs(subj,1,0,0);
        
        cd(params.savedir);
                
        % first check if we need to compute power for a subj/session
        doPow = 1;             
        for sess = 1:length(sessions)
            sessDir = fullfile(params.savedir,subj,[task '_events'],num2str(sessions(sess)));
            if ~exist(sessDir,'dir')
                doPow = 1;
                cd_mkdir(sessDir);
                
            else
                cd(sessDir);
                % if any electrode files are missing, make power for
                % session
                
                for iElec = 1:length(tal)
                    if strcmp(params.pow.type,'wavelet') || strcmp(params.pow.type,'hilbert')
                        if params.doBipol
                            fname = sprintf('%d-%d%s.mat',tal(iElec).channel(1),tal(iElec).channel(2),fileExt);
                        else
                            fname = sprintf('%d%s.mat',tal(iElec).channel(1),fileExt);
                        end
                    elseif strcmp(params.pow.type,'fft_slep')
                        fname = sprintf('%d-%d%s.mat',tal(iElec).channel(1),tal(iElec).channel(2),fileExt);                        
                    else
                        fprintf('Currently unsupported power type: %s.\n',params.pow.type);
                        return
                    end
                    if ~exist(fname,'file')
                        doPow = 1;
                        break;
                    end
                end
            end

            % compute power
            if doPow    
                
%                  if isempty(gcp('nocreate'))
%                      num_nodes = 30;mem = '12G';
%                      open_rhino2_pool(num_nodes,mem);
%                  end

                sess_events = events([events.session]==sessions(sess));
                parfor iElec = 1:length(tal)
                    ComputePow_local(tal(iElec),sess_events,sessDir,params,task,fileExt);
                end
                cd(sessDir);
                save('params.mat','params');                
            end
        end
        
    catch e
        cd_mkdir(fullfile(params.savedir,'errors'));
        dateFormat = 29; % output date as YYYY-MM-DD
        fname = sprintf('%s_%s_%s.mat',strrep(datestr(now,dateFormat),' ','_'), subj, task);
        save(fname,'e');
    end    
    display(sprintf('completed %s',subj));
    
end
% delete(gcp('nocreate'));




% ComputePow_local: parallel power over electrodes
function [] = ComputePow_local(Elec,events,sessDir,params,task,fileExt)

if params.doBipol
    for k=1:length(events)
        events(k).eegfile=regexprep(events(k).eegfile,'eeg.reref','eeg.noreref');
    end    
end

if strcmp(params.pow.type,'hilbert')
    [phase,amp,eeg] = gethilbertphase_bipol(Elec.channel,events,params.eeg.durationMS,params.eeg.offsetMS,params.eeg.bufferMS,params.pow.freqs,60);

    pow = amp.^2;
    pow(pow<=0) = eps;
    pow = log10(pow);    
   
    resampledrate = GetRateAndFormat(events(2));
    dsamp = round(1000/params.pow.timeStep);
    dmate = round(resampledrate/dsamp);
    
%     dsDur = ceil(size(pow,3)/params.pow.timeStep);    
    dsDur = ceil((size(phase,3) / (resampledrate/1000))/params.pow.timeStep);
    
    dpow = zeros(size(pow,1),size(pow,2),dsDur);   
    dphase = zeros(size(phase,1),size(phase,2),dsDur);
    
    for e = 1:size(amp,1)
        fprintf('%d ',e);
        for f = 1:size(amp,2)
            powTmp = decimate(double(pow(e,f,:)),dmate);
            dpow(e,f,1:length(powTmp)) = powTmp;
            phaseTmp = mod(decimate(double(unwrap(phase(e,f,:))),dmate)+pi,2*pi)-pi;            
            dphase(e,f,1:length(phaseTmp)) = phaseTmp;
        end
    end
    
    phase = permute(dphase,[2 3 1]);
    binPowMat = permute(dpow,[2 3 1]);

elseif ~params.useGetPhasePow
    phase = [];
    [EEG] = ComputeEEG(Elec.channel,events,params);
    [PowMat,~] = ComputePow(EEG,params);
    
    % fname = sprintf('%d-%d_raw.mat',Elec.channel(1),Elec.channel(2));
    cd_mkdir(sessDir);
    %save(fname,'PowMat');
    
    if strcmp(params.pow.type,'wavelet')
        [binPowMat] = BinPow(PowMat,params.pow.timeWin,params.pow.timeStep,...
            params.eeg.sampFreq,params.pow.freqs,params.pow.freqBins);
    elseif strcmp(params.pow.type,'fft_slep')
        binPowMat = mean(PowMat,2);
    end
    clear PowMat;
else
    dsamp = round(1000/params.pow.timeStep);
    [phase,binPowMat] = getphasepow_bipol(Elec.channel,events,params.eeg.durationMS,params.eeg.offsetMS,params.eeg.bufferMS,'freqs',params.pow.freqs,'width',params.pow.wavenum,'filtfreq',params.eeg.filtfreq,'filttype',params.eeg.filttype,'filtorder',params.eeg.filtorder,'downsample',dsamp);
    binPowMat = permute(binPowMat,[2 3 1]);
    binPowMat = log10(binPowMat);
    phase     = permute(phase,[2 3 1]);
end


% use mean/sd across word events to zscore
if strcmp(task,'RAM_YC1')
    baseInds    = params.eventsYC1(events);
elseif strcmp(task,'RAM_YC2')
    baseInds    = params.eventsYC2(events);
    if strcmp(params.pow.type,'fft_slep')
        baseInds    = params.eventsYC1(events);
    end
end
PowMean     = nanmean(nanmean(binPowMat(:,:,baseInds),3),2);
PowSTD      = std(nanmean(binPowMat(:,:,baseInds),2),[],3);
zPowMat     = zScorePow(binPowMat,PowMean,PowSTD);

sessOutput.pow = zPowMat;
sessOutput.meanBasePow = PowMean;
sessOutput.stdBasePow = PowSTD;
sessOutput.phase = phase;

if params.regressTrialNumber
    binPowRegress = NaN(size(binPowMat));
    x = [1:size(binPowMat,3)]';
    x = (x - mean(x))/std(x);
    for f = 1:size(binPowMat,1)
        for t = 1:size(binPowMat,2)
            y = permute(binPowMat(f,t,:),[3 1 2]);
            s = regstats(y,x,'linear',{'r'});
            binPowRegress(f,t,:) = s.r;
        end
    end
    
    
    PowMean     = nanmean(nanmean(binPowRegress(:,:,baseInds),3),2);
    PowSTD      = std(nanmean(binPowRegress(:,:,baseInds),2),[],3);
    zPowMat     = zScorePow(binPowRegress,PowMean,PowSTD);
    
    sessOutput.powCorr         = zPowMat;
    sessOutput.meanBasePowCorr = PowMean;
    sessOutput.stdBasePowCorr  = PowSTD;
end

if params.doBipol
    fname = sprintf('%d-%d%s.mat',Elec.channel(1),Elec.channel(2),fileExt);
else
    fname = sprintf('%d%s.mat',Elec.channel(1),fileExt);
end
cd_mkdir(sessDir);
save(fname,'sessOutput');


