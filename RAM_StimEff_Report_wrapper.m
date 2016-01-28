function [params,DidPow] = RAM_StimEff_Report_wrapper(Subject,Task,SessNum,powtype,ArtifactCorrect,ForceSess)
% function [params,DidPow] = RAM_StimEff_Report_wrapper(Subject,Task,SessNum,powtype,ArtifactCorrect,ForceSess)
%
% Wrapper function to do the stimulation classification analysis, and to
% create a .pdf report of the results.
%
% Last updated:
%       03/04/15    - YE added output for DidPow
%       02/26/15    - YE added section for collapsing individual session
%                   power across the same patient-stim-location (PSL)
%       02/13/15    - YE added input for powtype
%       01/21/15    - YE added check to skip if session has already been
%                   computed, with an optional input argument to force
%                   session to be re-computed
%       01/05/15    - YE added call to RAM_StimEff_RemoveArtifactPeaks, which
%                   implements the artifact removal procedures
%                   - YE updated so it reads for RAM_PatientInfo.xlsx to
%                   determine stimulation parameters that were used
%
%       12/07/14    - YE renamed to RAM_StimEff_Report_wrapper; also
%                   reorganized to call StimEff subfunctions and to
%                   operate on PAL/YC data

%% initialize base directory and parallel pool settings
basedir = GetBaseDir;

if ~exist('ForceSess','var')
    ForceSess = 0;
end

if isdir('/Users/yezzyat/Lab/data/events/') % running on bronx
    if isempty(gcp('nocreate'))
        parpool(10);
    end
elseif isdir('/data/events/') % running on rhino
    if matlabpool('SIZE')==0
        num_nodes = 15;mem = '4G';
        open_rhino_pool(num_nodes,mem);
    end
end

ClassParams = RAM_StimEff_Params(Subject,Task,powtype);

% load events, bipolar struct and indices for event types
events.all = get_sub_events(Task,Subject);
bpElec = getBipolarSubjElecs(Subject);
[inds] = RAM_GetWordAndStimInds(events.all,Task);
events.word = events.all(inds.word.pres);


%% compute power for all electrodes for this session
if ~exist('SessNum','var')
    SessList = unique([events.word.session]);
else
    SessList = SessNum;
end

% loop over sessions to compute power
for iSess = 1:length(SessList)
    if SessList(iSess) < 10
        SessSaveDir = fullfile(ClassParams.PowSaveDir,['Sess0',num2str(SessList(iSess))]);
    else
        SessSaveDir = fullfile(ClassParams.PowSaveDir,['Sess',num2str(SessList(iSess))]);
    end
    
    % check if session folder already created
    if exist(SessSaveDir,'dir')
        % check if all electrodes in tal struct have power
        DoSess = 0;
        cd(SessSaveDir);
        for iElec = 1:length(bpElec)
            if ~exist(sprintf('%d-%d_Pow_bin_zs.mat',bpElec(iElec).channel(1),bpElec(iElec).channel(2)),'file')
                DoSess = 1;
            end
        end
    else
        DoSess = 1;
    end
    %tic;
    if DoSess || ForceSess
        display(sprintf('COMPUTING POWER: %s\n',SessSaveDir));
        [params] = RAM_StimEff_ComputePow_PAR_MTaper(Subject,Task,SessList(iSess),SessSaveDir,powtype);
    else
        display(sprintf('Power found in %s\nSkipping session',SessSaveDir));
        cd(SessSaveDir);load params;
    end
    DidPow = DoSess || ForceSess;
    %toc;
end


%% nothing else to do in this function really, since you need info potentially collapsed across sessions.
%% --move the stat creation and report generation to a separate function.

%RAM_StimEff_CombinePSLPow(Subject,Task,params,powtype);














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% stuff below is from early version (~ Oct 2014) of stim classification code %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%
% % load the structure of session parameters. get the unique stim sessions for this task, and use to run classification within each parameter set individually
% SubjStimParams =  RAM_GetStimParamsFromDataBase(Subject);
% SessParams = RAM_DBStructToCell(SubjStimParams,Task);
% 
% % get the unique rows of the cell array of session parameters. the # of
% % unique rows will be the looping index for classification and
% % report-making.
% [UniqueSessParams,ndx,pos] = uniqueRowsCA(SessParams(:,2:5));
% 
% % set up the list of cross-validation units and indices of stim/no-stim
% % trials
% switch ClassParams.CrossVal
%     case 'Word'
%         % need to implement
%         %XValUnits = ?;
%         stim    = inds.word.pres & inds.word.stim & inds.word.PairPos==1;
%         nostim  = inds.word.pres & ~inds.word.stim & inds.word.PairPos==1;
%         
%     case 'List'
%         % hack for YC
% %         XValUnits = unique([[events.word.session]' [events.word.block]'],'rows');
% %         stim = inds.word.stim;
% %         nostim = inds.word.nostim;
%         XValUnits = unique([[events.word.session]' [events.word.list]'],'rows');
%         stim            = inds.word.pres & inds.word.stim & ~inds.list.nostim & inds.word.PairPos==1;
%         nostim          = inds.word.pres & inds.list.nostim & inds.word.PairPos==1;
%         nostimwords     = inds.list.stim & inds.word.nostim;
%         
%     case 'Session'
%         % need to implement
%         %XValUnits = ?;
% end
% 
% 
% %% Do classification, looping over unique parameter sets
% 
% % loop over unique parameter sets
% pcounter = 1;
% for iSet = 1:size(UniqueSessParams,1)
%     
%     if pcounter < 10
%         CurrClassSaveDir = [ClassParams.ClassSaveDir 'ParamSet0' num2str(pcounter)];
%     else
%         CurrClassSaveDir = [ClassParams.ClassSaveDir 'ParamSet' num2str(pcounter)];
%     end
%     
%     % if classification hasn't happened yet for this set of parameters
%     if exist(CurrClassSaveDir,'dir')
%         display(sprintf('** classification already performed for: Amp[%1.2f], Freq[%3.2f], Loc[%s]\n*** recomputing with all currently available data',...
%             UniqueSessParams{iSet,1},UniqueSessParams{iSet,2},UniqueSessParams{iSet,3}));
%     else
%         cd_mkdir(CurrClassSaveDir);
%     end
%     
%     % Determine the indices of the events to be included in this iteration
%     % of the loop
%     
%     % Get the session numbers for the sessions that matched the current
%     % parameter set. Figure out which x-validation units should be included
%     % in the current iteration.
%     
%     % hack for YC
% %    XValToKeep = logical(ones(length(XValUnits),1));
% 
%     switch ClassParams.CrossVal
%         case 'Word'
%             % need to implement
%         case 'List'
%             CurrParams = SessParams(pos==iSet,:);
%             XValToKeep = [];
%             TMP = zeros(size(XValUnits,1),1);
%             for iP = 1:size(CurrParams,1)
%                 if sum(isnan(CurrParams{iP,6}))
%                     TMP = TMP | XValUnits(:,1)==CurrParams{iP,1};
%                 else
%                     TMP = TMP | (XValUnits(:,1)==CurrParams{iP,1}) & ismember(XValUnits(:,2),CurrParams{iP,6});
%                 end
%             end
%             XValToKeep = logical(TMP);
%             
%         case 'Session'
%             % need to implement
%     end
%     
%     
%     % Add a check to make sure that the # of potential cross-validations is
%     % reasonable. For e.g. if only one list was stimulated at the current
%     % parameters, don't run the classification
%     if sum(XValToKeep) > 5
%         DoClass = 1;
%     else
%         DoClass = 0;
%     end
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%    classify the post-stim interval    %%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     if DoClass
%         
%         
%         
%         % get word onset bin and stim offset bin
%         clear PowMat
%         SessToIncl = unique([CurrParams{:,1}]);
%         [PowMat,params] = RAM_StimEff_LoadPow(bpElec(1),fullfile(ClassParams.PowSaveDir,['Sess0',num2str(SessToIncl(1))]));
%         %[PowMat,params] = RAM_StimEff_LoadPow(bpElec(1),ClassParams.PowSaveDir);
%         tBinStart = [params.eeg.offsetMS:params.pow.timeStep:params.eeg.durationMS-params.pow.timeWin];
%         tBinEnd = tBinStart+params.pow.timeWin;
%         StimOffTime = 4400; % time post-word onset at which stim turns off
%         StimOffBin = find((tBinStart<StimOffTime) & (tBinEnd>StimOffTime),1);
%         fInds = 1:size(PowMat,1);
%         tInds = StimOffBin:size(PowMat,2); 
%         
%         % run the classification on all electrodes
%         display(sprintf(' starting post-stim classification: param set %d of %d',iSet,size(UniqueSessParams,1)));
%         parfor iElec = 1:length(bpElec)
%         %for iElec = 1:length(bpElec)
%             %[~,~,~] = RAM_StimEff_SVM_SplitCrossVal(bpElec(iElec),events,inds,fInds,tInds,ClassParams.CrossVal,XValToKeep,Method,NoStimType,[],ClassParams.PowSaveDir,CurrClassSaveDir);
%             [~,~,~] = RAM_StimEff_SVM_SplitCrossVal(bpElec(iElec),events.word,stim,nostim,fInds,tInds,XValToKeep,ClassParams,CurrClassSaveDir);
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%    classify the mid-stim interval     %%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         % get word offset bin for first word (for FR2)
%         %[PowMat,~,params] = RAM_StimEff_LoadPow(bpElec(1),ClassParams.PowSaveDir);
%         tBinStart = [params.eeg.offsetMS:params.pow.timeStep:params.eeg.durationMS-params.pow.timeWin];
%         %tBinEnd = tBinStart+params.pow.timeWin;
%         Word1OffTime = 1600; % time post-word onset at which stim turns off
%         Word1OffBin = find((tBinStart>=Word1OffTime),1); % first bin after word offset
%         fInds = 1:size(PowMat,1);
%         tInds = Word1OffBin:Word1OffBin+99; clear PowMat;
%         
%         % run classification on all electrodes
%         display(sprintf(' starting mid-stim classification: param set %d of %d',iSet,size(UniqueSessParams,1)));
%         parfor iElec = 1:length(bpElec)
%             %[~,~,~] = RAM_StimEff_SVM_SplitCrossVal(bpElec(iElec),events,inds,fInds,tInds,CrossVal,XValToKeep,Method,NoStimType,[],ClassParams.PowSaveDir,CurrClassSaveDir);
%             [~,~,~] = RAM_StimEff_SVM_SplitCrossVal(bpElec(iElec),events.word,stim,nostim,fInds,tInds,XValToKeep,ClassParams,CurrClassSaveDir);
%         end
%         
%         
%         
%         % save a record of the parameters and sessions that went into this
%         % iteration of classification
%         cd(CurrClassSaveDir);
%         ParamsSessData = SessParams(pos==iSet,:);
%         save ParamsSessData;
%         
%         pcounter = pcounter+1; % keep track of how many parameter sets were actually used for classification
%         
%     end
%     
%    
% end % iSet -- unique stimulation parameter sets
% 
% 
% %% Get mean classifier performance for mid-stim and post-stim windows
% SessToIncl = unique([CurrParams{:,1}]);
% [PowMat,params] = RAM_StimEff_LoadPow(bpElec(1),fullfile(ClassParams.PowSaveDir,['Sess0',num2str(SessToIncl(1))]));
% %for iSet = 1:size(UniqueSessParams,1)
% for iSet = 1:2
%     
%     %if strcmp(Subject,'R1003P') && strcmp(Task,'RAM_PAL2')
%         if iSet < 10
%             CurrClassSaveDir = [ClassParams.ClassSaveDir 'ParamSet0' num2str(iSet)];
%         else
%             CurrClassSaveDir = [ClassParams.ClassSaveDir 'ParamSet' num2str(iSet)];
%         end
%     %else
%     %    CurrClassSaveDir = ClassParams.ClassSaveDir;
%     %end
%         
%         if isdir(CurrClassSaveDir)
%             
%             % get mean classifier performance across cross-validations
%             tBinStart = [params.eeg.offsetMS:params.pow.timeStep:params.eeg.durationMS-params.pow.timeWin];
%             tBinEnd = tBinStart+params.pow.timeWin;
%             StimOffTime = 4400; % time post-word onset at which stim turns off
%             StimOffBin = find((tBinStart<StimOffTime) & (tBinEnd>StimOffTime),1);
%             tInds = StimOffBin:size(PowMat,2);
%             
%             PostClassWeights = nan(length(bpElec),size(ClassParams.FreqBins,1),2);
%             nBins = 10; % number of bins over which to average to compute mean classifier performance
%             for iElec = 1:length(bpElec)
%                 cd(CurrClassSaveDir);
%                 load(sprintf('%d-%d_SVM_CV-%s_f2-200_t%d-%d.mat', bpElec(iElec).channel(1),bpElec(iElec).channel(2),ClassParams.CrossVal,tBinStart(tInds(1)),tBinEnd(tInds(end))));
%                 
%                 if iElec==1 & iSet==1
%                     PostStimPerfMat = nan(length(bpElec),length(Perf),size(UniqueSessParams,1));
%                 end
%                 PostStimPerfMat(iElec,:,iSet) = Perf;
%                 
%                 % get mean classifier weights for this electrode
%                 % weight matrix: nFreqs x n crossval x nTimeBins
%                 for iBin = 1:size(ClassParams.FreqBins,1)
%                     CurrFInds = params.pow.freqBins>=ClassParams.FreqBins(iBin,1) & params.pow.freqBins<ClassParams.FreqBins(iBin,2);
%                     PostClassWeights(iElec,iBin,iSet) = nanmean(nanmean(nanmean((W(CurrFInds,:,nBins)))));
%                 end
%                 
%             end
%             
%             % get mean classifier performance and feature weights across cross-validations
%             tBinStart = [params.eeg.offsetMS:params.pow.timeStep:params.eeg.durationMS-params.pow.timeWin];
%             %tBinEnd = tBinStart+params.pow.timeWin;
%             Word1OffTime = 1600; % time post-word onset at which stim turns off
%             Word1OffBin = find((tBinStart>=Word1OffTime),1); % first bin after word offset
%             fInds = 1:size(PowMat,1);
%             tInds = Word1OffBin:Word1OffBin+99;
%             
%             MidClassWeights = nan(length(bpElec),size(ClassParams.FreqBins,1),2);
%             for iElec = 1:length(bpElec)
%                 cd(CurrClassSaveDir);
%                 load(sprintf('%d-%d_SVM_CV-%s_f2-200_t%d-%d.mat', bpElec(iElec).channel(1),bpElec(iElec).channel(2),ClassParams.CrossVal,tBinStart(tInds(1)),tBinEnd(tInds(end))));
%                 if iElec==1 && iSet==1
%                     MidStimPerfMat = nan(length(bpElec),length(Perf),size(UniqueSessParams,1));
%                 end
%                 MidStimPerfMat(iElec,:,iSet) = Perf;
%                 
%                 % get mean classifier weights for this electrode
%                 % weight matrix: nFreqs x n crossval x nTimeBins
%                 for iBin = 1:size(ClassParams.FreqBins,1)
%                     CurrFInds = params.pow.freqBins>=ClassParams.FreqBins(iBin,1) & params.pow.freqBins<ClassParams.FreqBins(iBin,2);
%                     MidClassWeights(iElec,iBin,iSet) = nanmean(nanmean(nanmean((W(CurrFInds,:,:)))));
%                 end
%                 
%             end
%             
%         end % isdir(CurrClassSaveDir)
%     
% end % iSet

%%


%% get the inter-electrode distances (for use later)
% get indices of electrode groups
% %for iSet = 1:size(UniqueSessParams,1)
%     SubjStimParams =  RAM_GetStimParamsFromDataBase(Subject);
%     StimElecTags = unique({SubjStimParams.(Task).StimTag});
%     StimElecGrps = unique({SubjStimParams.(Task).StimGrp});
%     isDepth = inStruct(bpElec,'strcmp(eType,''D'')');
%     
%     % hack for YC
%     StimElecInd = inStruct(bpElec,'strcmp(tagName,''DA1-DA2'')');
%     %StimElecInd = inStruct(bpElec,sprintf('strcmp(tagName,''%s'')',StimElecTags{1}));
%     StimElecTalXYZ = [bpElec(StimElecInd).x bpElec(StimElecInd).y bpElec(StimElecInd).z];
%     StimElecIndXYZ = [bpElec(StimElecInd).indivSurf.x bpElec(StimElecInd).indivSurf.y bpElec(StimElecInd).indivSurf.z];
%     
%     AllElecTalXYZ = nan(length(bpElec),3);
%     AllElecIndXYZ = nan(length(bpElec),3);
%     for iElec = 1:length(bpElec)
%         % get spatial coordinates of each electrode
%         AllElecTalXYZ(iElec,:) = [bpElec(iElec).x bpElec(iElec).y bpElec(iElec).z];
%         AllElecIndXYZ(iElec,:) = [bpElec(iElec).indivSurf.x bpElec(iElec).indivSurf.y bpElec(iElec).indivSurf.z];
%     end
%     
%     % load SVM classifier performance at each electrode
%     TalDist = pdist([StimElecTalXYZ;AllElecTalXYZ]); TalDist = TalDist(1:length(bpElec));
%     IndDist = pdist([StimElecIndXYZ;AllElecIndXYZ]); IndDist = IndDist(1:length(bpElec));
% %end


%% load power and average for post-stim and mid-stim intervals. do separately for each parameter set
% % cd(ClassParams.ClassSaveDir);
% % d = dir('ParamSet*');
% % for iD = 2:length(d)
% %     cd(d(iD).name);
% %     load ParamsSessData
% %     TMP = zeros(length(events.word),1);
% %     for iP = 1:size(ParamsSessData,1)
% %         if sum(isnan(ParamsSessData{iP,6}))
% %             TMP = TMP | ([events.word.session]==ParamsSessData{iP,1})';
% %         else
% %             TMP = TMP | ([events.word.session]==ParamsSessData{iP,1} & ismember([events.word.list],ParamsSessData{iP,6}))';
% %         end
% %     end
% %     SetEvInds = TMP; clear TMP;
% SetEvInds = logical(ones(size(nostim')));
% %ClassParams.PowSaveDir = [ClassParams.PowSaveDir 'ParamSet02/'];
%     
%     PostStimPowMat = nan(length(bpElec),size(ClassParams.FreqBins,1),3);
%     MidStimPowMat = nan(length(bpElec),size(ClassParams.FreqBins,1),3);
%         DepthInds = find(isDepth);
%     for pElec = 1:length(DepthInds)
%         iElec = DepthInds(pElec);
%         
%         [PowMat,params] = RAM_StimEff_LoadPow(bpElec(iElec),ClassParams.PowSaveDir);
%         
%         % mid-stim
%         tInds = Word1OffBin:Word1OffBin+99;
%         for iBin = 1:size(ClassParams.FreqBins,1)
%             CurrFInds = params.pow.freqBins>=ClassParams.FreqBins(iBin,1) & params.pow.freqBins<ClassParams.FreqBins(iBin,2);
%             MidStimPowMat(iElec,iBin,1) = nanmean(nanmean(nanmean((PowMat(CurrFInds,tInds(1):tInds(end),stim' & SetEvInds)))));
%             %MidStimPowMat(iElec,iBin,2) = nanmean(nanmean(nanmean((PowMat(CurrFInds,tInds(1):tInds(end),nostimwords' & SetEvInds)))));
%             MidStimPowMat(iElec,iBin,2) = nanmean(nanmean(nanmean((PowMat(CurrFInds,tInds(1):tInds(end),nostim' & SetEvInds)))));
%         end
%         
%         % post-stim
%         %tInds = StimOffBin:size(PowMat,2);
%         tInds = StimOffBin:StimOffBin+25;
%         for iBin = 1:size(ClassParams.FreqBins,1)
%             CurrFInds = params.pow.freqBins>=ClassParams.FreqBins(iBin,1) & params.pow.freqBins<ClassParams.FreqBins(iBin,2);
%             PostStimPowMat(iElec,iBin,1) = nanmean(nanmean(nanmean((PowMat(CurrFInds,tInds(1):tInds(end),stim' & SetEvInds)))));
%             %PostStimPowMat(iElec,iBin,2) = nanmean(nanmean(nanmean((PowMat(CurrFInds,tInds(1):tInds(end),nostimwords' & SetEvInds)))));
%             PostStimPowMat(iElec,iBin,2) = nanmean(nanmean(nanmean((PowMat(CurrFInds,tInds(1):tInds(end),nostim' & SetEvInds)))));
%         end
%         clear PowMat
%     end
%     
%     % create the plots of stim-nostim power averaged within frequency band
% %     if ~isdir([ClassParams.PowPlotDir d(iD).name])
% %         mkdir([ClassParams.PowPlotDir d(iD).name]);
% %     end
%     % create plot for during stim
% 
%     %iDepthE = iElec;
%     for iDepthE = 1:length(DepthInds)
%         figure(1);clf;
%         bardata = [MidStimPowMat(DepthInds(iDepthE),:,1) MidStimPowMat(DepthInds(iDepthE),:,2) MidStimPowMat(DepthInds(iDepthE),:,3)];
%         barcolors = [repmat([0 0 0],size(ClassParams.FreqBins,1),1); repmat([.5 .5 .5],size(ClassParams.FreqBins,1),1);repmat([1 1 1],size(ClassParams.FreqBins,1),1)];
%         barlocs =  [1:4:size(ClassParams.FreqBins,1)*4 (1:4:size(ClassParams.FreqBins,1)*4)+1 (1:4:size(ClassParams.FreqBins,1)*4)+2];
%         for iBar = 1:numel(bardata)
%             hbar = bar(barlocs(iBar),bardata(iBar)); hold on;
%             set(get(hbar,'Children'),'FaceColor',barcolors(iBar,:));
%             set(get(hbar,'Children'),'LineWidth',2);
%         end
%         %hEB = errorbar(barlocs,[CatIRT BtwCatIRT],[nan nan],'Linestyle','none');
%         %set(hEB,'LineWidth',2); set(hEB,'Color',[0 0 0]); hold on;
%         
%         set(gca,'FontSize',14);
%         set(gca,'Xlim',[0 ceil(max(barlocs))+1])
%         set(gca,'ylim',[min(bardata)-(.1)*range(bardata) max(bardata)+(.3)*range(bardata)]);
%         xtickinds = [2:4:size(ClassParams.FreqBins,1)*4];
%         %xtickinds = [1.5:3:size(ClassParams.FreqBins,1)*3];
%         set(gca,'xtick',xtickinds);
%         format_ticks(gca,ClassParams.FreqBinLabels);
%         %[hleg hlegobj hlegout ~] = legend('Stim','No Stim Words','No Stim Lists','Location','NorthWest');
%         [hleg hlegobj hlegout ~] = legend('Stim','No Stim','Location','NorthWest');
%         set(hlegout(1),'LineWidth',2); set(hlegout(2),'LineWidth',2);
%         set(get(hlegobj(3),'Children'),'FaceColor',[0 0 0]);
%         set(get(hlegobj(4),'Children'),'FaceColor',[.5 .5 .5]);
%         %set(get(hlegobj(6),'Children'),'FaceColor',[1 1 1]);
%         set(hleg,'Box','off');
%         title(bpElec(DepthInds(iDepthE)).tagName);
%         box off
%         ylabel('Power (z)');
%         set(gcf,'PaperPositionMode','auto');
%         Figs.PowMid = sprintf('%s_Chan%d-%d_Pow_MidStim.eps',bpElec(DepthInds(iDepthE)).tagName,bpElec(DepthInds(iDepthE)).channel(1),bpElec(DepthInds(iDepthE)).channel(2));
%         %print(fullfile([ClassParams.PowPlotDir d(iD).name],Figs.PowMid));
%         print(fullfile([ClassParams.PowPlotDir],Figs.PowMid));
%         
%     end
%     
%     % create plot for post stim
%     for iDepthE = 1:length(DepthInds)
%         figure(1);clf;
%         bardata = [PostStimPowMat(DepthInds(iDepthE),:,1) PostStimPowMat(DepthInds(iDepthE),:,2) PostStimPowMat(DepthInds(iDepthE),:,3)];
%         barcolors = [repmat([0 0 0],size(ClassParams.FreqBins,1),1); repmat([.5 .5 .5],size(ClassParams.FreqBins,1),1);repmat([1 1 1],size(ClassParams.FreqBins,1),1)];
%         barlocs =  [1:4:size(ClassParams.FreqBins,1)*4 (1:4:size(ClassParams.FreqBins,1)*4)+1 (1:4:size(ClassParams.FreqBins,1)*4)+2];
%         for iBar = 1:numel(bardata)
%             hbar = bar(barlocs(iBar),bardata(iBar)); hold on;
%             set(get(hbar,'Children'),'FaceColor',barcolors(iBar,:));
%             set(get(hbar,'Children'),'LineWidth',2);
%         end
%         %hEB = errorbar(barlocs,[CatIRT BtwCatIRT],[nan nan],'Linestyle','none');
%         %set(hEB,'LineWidth',2); set(hEB,'Color',[0 0 0]); hold on;
%         
%         set(gca,'FontSize',14);
%         set(gca,'Xlim',[0 ceil(max(barlocs))+1])
%         set(gca,'ylim',[min(bardata)-(.1)*range(bardata) max(bardata)+(.3)*range(bardata)]);
%         xtickinds = [2:4:size(ClassParams.FreqBins,1)*4];
%         %xtickinds = [1.5:3:size(ClassParams.FreqBins,1)*3];
%         set(gca,'xtick',xtickinds);
%         format_ticks(gca,ClassParams.FreqBinLabels);
%         %[hleg hlegobj hlegout ~] = legend('Stim','No Stim Words','No Stim Lists','Location','NorthWest');
%         [hleg hlegobj hlegout ~] = legend('Stim','No Stim','Location','NorthWest');
%         set(hlegout(1),'LineWidth',2); set(hlegout(2),'LineWidth',2);
%         set(get(hlegobj(3),'Children'),'FaceColor',[0 0 0]);
%         set(get(hlegobj(4),'Children'),'FaceColor',[.5 .5 .5]);
%         %set(get(hlegobj(6),'Children'),'FaceColor',[1 1 1]);
%         set(hleg,'Box','off');
%         title(bpElec(DepthInds(iDepthE)).tagName);
%         box off
%         ylabel('Power (z)');
%         set(gcf,'PaperPositionMode','auto');
%         Figs.PowMid = sprintf('%s_Chan%d-%d_Pow_PostStim.eps',bpElec(DepthInds(iDepthE)).tagName,bpElec(DepthInds(iDepthE)).channel(1),bpElec(DepthInds(iDepthE)).channel(2));
%         %print(fullfile([ClassParams.PowPlotDir d(iD).name],Figs.PowMid));
%         print(fullfile([ClassParams.PowPlotDir],Figs.PowMid));
%         
%     end
%     
%     
% %     cd([ClassParams.PowPlotDir d(iD).name]);cd ../
% %     %unix(['chgrp -R RAM_eeg ' ClassParams.PowPlotDir]);
% %     
% % end % for iD
% 
% 
% %% write the latex file
% 
% % write the preamble
% fid = fopen(fullfile(RepSaveDir,ReportName),'w');
% if fid==-1;
%     error(sprintf('cannot open %s',ReportName))
% end
% 
% % write out preamble
% fprintf(fid,'\\documentclass[a4paper]{article} \n');
% fprintf(fid,'\\usepackage{graphicx,multirow} \n');
% fprintf(fid,'\\usepackage{epstopdf} \n');
% fprintf(fid,'\\usepackage{subfigure,amsmath} \n');
% fprintf(fid,'\\usepackage{wrapfig} \n');
% fprintf(fid,'\\usepackage{pdfpages}\n');
% fprintf(fid,'\\usepackage{mathtools}\n');
% fprintf(fid,'\\usepackage{array}\n');
% %fprintf(fid,'\\usepackage[table]{xcolor} \n');
% fprintf(fid,'\n');
% fprintf(fid,'\\addtolength{\\oddsidemargin}{-.875in} \n');
% fprintf(fid,'\\addtolength{\\evensidemargin}{-.875in} \n');
% fprintf(fid,'\\addtolength{\\textwidth}{1.75in} \n');
% fprintf(fid,'\\addtolength{\\topmargin}{-.75in} \n');
% fprintf(fid,'\\addtolength{\\textheight}{1.75in} \n');
% fprintf(fid,'\n');
% fprintf(fid,'\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}} \n');
% 
% % Start the document
% fprintf(fid,'\\begin{document}\n\n\n');
% 
% % write start of the report tex file
% [fid, msg] = fopen(fullfile(RepSaveDir,ReportName),'a');
% if fid==-1;
%     error(sprintf('cannot open %s: %s',ReportName,msg))
% end
% 
% %blurb_TexName = 'RAM_FR_blurb.tex';
% fprintf(fid,'\\section*{%s %s Stimulation Classification Summary}\n\n',regexprep(Subject,'_','\\_'),regexprep(Task,'_',' '));
% fprintf(fid,'\\vspace{0.3 in}\n\n');
% 
% % create the depths table for data DURING stimulation
% %   elec pair / class % during stim / class % post stim / mean classifier feature weights
% fprintf(fid,['\\begin{table}[!ht] \n',...
%     '\\centering \n',...
%     '\\begin{tabular}{C{3cm} | C{2 cm} | C{1.25cm} C{1.25cm} C{1.25cm} C{1.25cm} C{1.25cm} C{1.25cm}}\n' ...
%     '\\multicolumn{8}{c}{\\textbf{\\large During Stimulation}} \\\\ \n' ...
%     '\\hline \n',...
%     'Electrode Pair & Mean Hit & $\\theta$ & $\\alpha$ & $\\beta$ & low $\\gamma$ & high $\\gamma$ & HFA \\\\ \n',...
%     '\\hline \n ']);
% DepthTableData = [mean(MidStimPerfMat(isDepth,:),2),...
%     MidClassWeights(isDepth,:)];
% for iElec = 1:length(DepthInds)
%     fprintf(fid,'%s & %0.2f & %2.1f & %2.1f & %2.1f & %2.1f & %2.1f & %2.1f \\\\ \n', bpElec(DepthInds(iElec)).tagName,...
%         DepthTableData(iElec,1),...
%         DepthTableData(iElec,2),DepthTableData(iElec,3),...
%         DepthTableData(iElec,4),DepthTableData(iElec,5),...
%         DepthTableData(iElec,6),DepthTableData(iElec,7));
% end
% fprintf(fid,'\\hline \n');
% fprintf(fid,'\\end{tabular}\n');
% fprintf(fid,['\\caption{',...
%     'Classifier data shown in depth electrodes during stimulation. Col 2: mean performance for the ISI between stimulated words. ',...
%     'Col 3-8: mean classifier weight for features in each frequency band. Positive weights indicate features were predictive of stimulation. ',...
%     'Negative weights indicate features were predictive of non-stimulation.',...
%     '$\\theta$ = %d-%d Hz; $\\alpha$ = %d-%d Hz; $\\beta$ = %d-%d Hz; low $\\gamma$ = %d-%d Hz; high $\\gamma$ = %d-%d Hz; HFA = %d-%d Hz. }\n'],...
%     FreqBins(1,1),FreqBins(1,2),...
%     FreqBins(2,1),FreqBins(2,2),...
%     FreqBins(3,1),FreqBins(3,2),...
%     FreqBins(4,1),FreqBins(4,2),...
%     FreqBins(5,1),FreqBins(5,2),...
%     FreqBins(6,1),FreqBins(6,2));
% 
% fprintf(fid,'\\end{table}\n\n\n');
% 
% 
% % create the depths table for data POST stimulation 
% %   elec pair / class % during stim / class % post stim / mean classifier feature weights
% fprintf(fid,['\\begin{table}[!ht] \n',...
%         '\\centering \n',...
%         '\\begin{tabular}{C{3cm} | C{2cm} | C{1.25cm} C{1.25cm} C{1.25cm} C{1.25cm} C{1.25cm} C{1.25cm}}\n' ...
%         '\\multicolumn{8}{c}{\\textbf{\\large Post Stimulation}} \\\\ \n' ...
%         '\\hline \n',...
%         'Electrode Pair & Mean Hit & $\\theta$ & $\\alpha$ & $\\beta$ & low $\\gamma$ & high $\\gamma$ & HFA \\\\ \n',...
%         '\\hline \n ']);
% DepthTableData = [mean(PostStimPerfMat(isDepth,:),2),...
%     PostClassWeights(isDepth,:)];
% for iElec = 1:length(DepthInds)
%     fprintf(fid,'%s & %0.2f & %2.1f & %2.1f & %2.1f & %2.1f & %2.1f & %2.1f \\\\ \n', bpElec(DepthInds(iElec)).tagName,...
%         DepthTableData(iElec,1),...
%         DepthTableData(iElec,2),DepthTableData(iElec,3),...
%         DepthTableData(iElec,4),DepthTableData(iElec,5),...
%         DepthTableData(iElec,6),DepthTableData(iElec,7));
% end
% fprintf(fid,'\\hline \n');
% fprintf(fid,'\\end{tabular}\n');
% fprintf(fid,['\\caption{',...
%     'Classifier data shown in depth electrodes following stimulation offset. Col 2: mean performance for 0-100ms of the ISI following stimulation offset. ',...
%     'Col 3-8: mean classifier weight for features in each frequency band. Positive weights indicate features were predictive of stimulation. ',...
%     'Negative weights indicate features were predictive of non-stimulation. Frequency bands same as above.}\n']);
% fprintf(fid,'\\end{table}\n\n\n');
% 
% fprintf(fid,'\\clearpage\n\n');
% 
% 
% % write out the power plots
% fprintf(fid,[...
%   '\\begin{figure}[t] \n' ... 
%   '\\centering \n' ...
%   '\\begin{tabular}{cc} \n' ...
%   '\\multicolumn{2}{c}{\\textbf{\\large Mean Power}} \\vspace{0.2 in} \\\\ \n',...
%   '\\vspace{0.1 in}\n',...
%    'During Stimulation & Post Stimulation \\\\ \n']);
% 
% FigScale = .5;
% for iElec = 1:length(DepthInds)
%     FigMidName = fullfile(PowPlotDir,sprintf('%s_Chan%d-%d_Pow_MidStim.eps',bpElec(DepthInds(iElec)).tagName,bpElec(DepthInds(iElec)).channel(1),bpElec(DepthInds(iElec)).channel(2)));
%     FigPostName = fullfile(PowPlotDir,sprintf('%s_Chan%d-%d_Pow_PostStim.eps',bpElec(DepthInds(iElec)).tagName,bpElec(DepthInds(iElec)).channel(1),bpElec(DepthInds(iElec)).channel(2)));
%     fprintf(fid,'\\includegraphics[scale=%0.2f]{%s} & \\includegraphics[scale=%0.2f]{%s} \\\\ \n',FigScale,FigMidName,FigScale,FigPostName);
%     fprintf(fid,'\\hline \n\n');
%     
%     % if up to 4th row, start a new page
%     if (mod(iElec,4)==0) & (iElec < length(DepthInds))
%         fprintf(fid,['\\end{tabular}\n',...
%             '\\end{figure} \n\n']);
%         fprintf(fid,'\\clearpage\n\n');
%         fprintf(fid,[...
%             '\\begin{figure}[t] \n' ... 
%             '\\centering \n' ...
%             '\\begin{tabular}{cc} \n' ...
%             '\\multicolumn{2}{c}{\\textbf{\\large Mean Power}} \\vspace{0.2 in} \\\\ \n',...
%             '\\vspace{0.1 in}\n',...
%             'During Stimulation & Post Stimulation \\\\ \n']);
%     end
% end
% 
% fprintf(fid,['\\end{tabular}\n',...
%     '\\end{figure} \n\n']);
% 
% % write end of the report tex file
% fid = fopen(fullfile(RepSaveDir,ReportName),'a');
% if fid==-1;
%     error(sprintf('cannot open %s',ReportName))
% end
% fprintf(fid,'\\end{document}\n\n\n');
% fclose(fid);
% 
% % compile to pdf
% curr_dir = pwd;
% cd(RepSaveDir);
% unix(['pdflatex -shell-escape ' fullfile(RepSaveDir, ReportName)]);
% unix(['rm ' ReportName(1:end-3) 'aux']);
% unix(['rm ' ReportName(1:end-3) 'log']);
% cd (RepSaveDir);
% unix(['chgrp RAM_eeg ' ReportName(1:end-3) '*']);
% 


%% get mean classifier weights during stim and post-stim



%% get stim/no stim power difference in all freq bands



%% get # peaks in each freq band for stim/nostim


