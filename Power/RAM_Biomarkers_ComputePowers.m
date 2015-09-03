function [] = RAM_Biomarkers_ComputePowers(Task,params)
% function [] = RAM_Biomarkers_ComputePowers(Task,params)
%
% Function to compute powers for multivariate biomarker analysis.
% Supersedes the use of RAM_reporting-generated powers.
%
% Inputs:
%       Task        'RAM_FR1'
%
%       params      params struct created from CreateEEGAndPowParams. Must
%                   specify params.events as a function handle:
%                       e.g. params.events = @(events)strcmp({events.type},'WORD')
%                   Also, must specify a location to check/put saved
%                   powers:
%                       e.g. '/data10/scratch/yezzyat/RAM/biomarker/power/'
%
% Last updated:
%       07/23/15    Created function


SubjList = get_subs(Task);
KeepSs = strncmp(SubjList,'R1',2) | strncmp(SubjList,'TJ',2) | strncmp(SubjList,'UT',2);
SubjList = SubjList(KeepSs);

for iSubj = 1:length(SubjList)
    
    try
        doPow = 0;
        
        Subject = SubjList{iSubj};
        events = get_sub_events(Task,Subject);
        %events = params.events(events);
        SessList = unique([events.session]);
        bpStruct = getBipolarSubjElecs(Subject,1,1,0);
        
        cd(params.savedir);
        % check for session power folders.
        for iSess = 1:length(SessList)
            SessPowDir = fullfile(params.savedir,Subject,[Task '_events'],num2str(SessList(iSess)));
            if ~exist(SessPowDir,'dir')
                doPow = 1;
                cd_mkdir(SessPowDir);
            else
                cd(SessPowDir);
                % check individual electrode files
                for iElec = 1:length(bpStruct)
                    fname = sprintf('%d-%d.mat',bpStruct(iElec).channel(1),bpStruct(iElec).channel(2));
                    if ~exist(fname,'file')
                        doPow = 1;
                        break;
                    end
                end
            end
            
            % compute power
            if doPow
                [~,hostname] = system('hostname');
                if ~isempty(strfind(hostname,'bronx.psych.upenn.edu'))
                    if isempty(gcp('nocreate'))
                        parpool(10);
                    end
                else
                    %elseif ~isempty(strfind(hostname,'rhino2'))
                    if isempty(gcp('nocreate'))
                        num_nodes = 30;mem = '12G';
                        open_rhino2_pool(num_nodes,mem);
                    end
                end
                
                SessEvents = events([events.session]==SessList(iSess));
                parfor iElec = 1:length(bpStruct)
                    ComputePow_local(bpStruct(iElec),SessEvents,SessPowDir,params);
                end
                cd(SessPowDir);
                save('params.mat','params');
                
            end
        end
        
    catch e
        cd_mkdir(fullfile(params.savedir,'errors'));
        dateFormat = 29; % output date as YYYY-MM-DD
        fname = sprintf('%s_%s_%s.mat',strrep(datestr(now,dateFormat),' ','_'), Subject, Task);
        save(fname,'e');
    end
    
    display(sprintf('completed %s',Subject));
    
end % iSubj

delete(gcp('nocreate'));




%% ComputePow_local: parallel power over electrodes
function [] = ComputePow_local(Elec,events,SessPowDir,params)

[EEG] = ComputeEEG(Elec.channel,events,params);
[PowMat,~] = ComputePow(EEG,params);
fname = sprintf('%d-%d_raw.mat',Elec.channel(1),Elec.channel(2));
cd_mkdir(SessPowDir);
%save(fname,'PowMat');

[binPowMat] = BinPow(PowMat,params.pow.timeWin,params.pow.timeStep,...
    params.eeg.sampFreq,params.pow.freqs,params.pow.freqBins);

clear PowMat;

% use mean/sd across word events to zscore
baseInds    = params.events(events);
PowMean     = nanmean(nanmean(binPowMat(:,:,baseInds),3),2);
PowSTD      = std(nanmean(binPowMat(:,:,baseInds),2),[],3);

zPowMat = zScorePow(binPowMat,PowMean,PowSTD);

sessOutput.pow = zPowMat;
sessOutput.meanBasePow = PowMean;
sessOutput.stdBasePow = PowSTD;

fname = sprintf('%d-%d.mat',Elec.channel(1),Elec.channel(2));
cd_mkdir(SessPowDir);
save(fname,'sessOutput');


