function YC1_weightsByRegion(subjs,params,overwrite)
% function YC1_weightsByRegion(subjs,params,overwrite)
%
% Inputs:
%
%      subjs - cell array of subject strings (default get_subs('RAM_YC1'))
%     params - params structure (default is returned by multiParams)
%  overwrite - boolean. if true, overwrite existing figures
%
% Make a report of the classifier weights for each subject in YC1 that has
% classification and chance clasification data. Also make a group average
% report.
%
% For each subject, figures are:
%
%   - Classifier weights for each electrode and frequency, averaged over
%     time
%   - Classifier weights for each electrode and time, averaged over
%     frequencies
%   - Classifier weights for for time and frequency, averaged over
%     electrodes
%   - A series of figures with the weights broken down by specific brain
%     region
%
%
% For the group report, there is time x frequency spectrogram, as well as
% averaged over time and over frequency. Also includes brain region average
% weights
%
%  Reports are saved to <location of subject
%  data>/reports/lassoWeights_report.pdf and group_lassoWeights_report.pdf

% if not given, use default params
if ~exist('params','var') || isempty(params)
    params = multiParams();
end


% overwrite exisiting figures?
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = true;
end

% save directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
prefix = '';
if ismac
    prefix = '/Volumes/data';
end

dataDir = fullfile(prefix,params.basePath,f(params.useCorrectedPower,y));
saveDir = fullfile(dataDir,'reports');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% figure directory
figDir = fullfile(saveDir,'figs');
if ~exist('figDir','dir')
    mkdir(figDir)
end

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end
% subjs = subjs(~strcmp(subjs,'R1061T'))

nTimes = size(params.timeBins,1);
nFreqs = size(params.freqBins,1);
spectTimeFreq             = NaN(nFreqs,nTimes,length(subjs));
spectTimeFreqNonAbs       = NaN(nFreqs,nTimes,length(subjs));


figs = [];
for s = 1:length(subjs)
    subj = subjs{s};
    fprintf('Processing %s.\n',subjs{s})
    
    figs_subj = struct('subj',[],'Spect_Freq',[],'Spect_Time',[],...
        'Spect_TxF',[],'Region_Bar',[],'nElecs',[]);
    
    % load the weights for this subject
    out = loadWeightsByRegion(subjs{s},dataDir);
    if ~isempty(out)
        
        if ~exist('regionDataAll','var')
            nRegions                   = length(out.regionsAll);
            regionDataAll              = NaN(nRegions,length(subjs),nTimes);
            regionDataNonZeroAll       = NaN(nRegions,length(subjs),nTimes);
            regionDataAllFreqsTimes    = NaN(nRegions,length(subjs),nTimes,nFreqs);
            regionDataAllFreqsBestTime = NaN(nRegions,length(subjs),nFreqs);
        end
        
        
        figs_subj.subj = subj;
        figs_subj.Region_Bar = {};        
        
        % for each subject, we will plot:
        %    elec x freq spectrogram
        %    elec x time spectrogram
        %    time x freq spectrogram
        
        %% FIGURE 1 - elec x freq (averaged across time)
        fname = fullfile(figDir,[subj 'Spect_Freq.png']);
        figs_subj.Spect_Freq = fname;
        
        % create y labels based on freq bins
        yBinsStrF = {};
        for x = 1:size(params.freqBins,1)
            yBinsStrF{x} = [num2str(params.freqBins(x,1)), '-', num2str(params.freqBins(x,2))];
        end
        
        if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
            figure(1)
            clf
            imagesc(squeeze(mean(abs(out.meanWeightsPerTimeSort),3)))
            colormap('jet')
            colorbar
            set(gca,'ytick',1:nFreqs)
            set(gca,'yticklabel',yBinsStrF)
            set(gca,'xtick',[])
            set(gca,'fontsize',16)
            set(gca,'xtick',out.regionCutoffs);
            set(gca,'xticklabel','|');
            set(gca,'xticklabel',out.regions);
            set(gca,'XTickLabelRotation',270)
            xlabel('Electrodes','fontsize',20)
            ylabel('Frequency','fontsize',20)
            print('-dpng','-loose',fname);
        end
        
        
        %% FIGURE 2 - elec x time (averaged across freq)
        fname = fullfile(figDir,[subj 'Spect_Time.png']);
        figs_subj.Spect_Time = fname;
        
        % create y labels based on freq bins
        yBins    = round((params.timeBins) / 10) / 100;
        yBinsStrT = {};
        for x = 1:size(yBins,1)
            yBinsStrT{x} = [num2str(yBins(x,1)), '-', num2str(yBins(x,2))];
        end
        
        if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
            figure(2)
            clf
            imagesc(squeeze(mean(abs(out.meanWeightsPerTimeSort),1))');
            colormap('jet')
            colorbar
            set(gca,'ytick',1:nTimes)
            set(gca,'yticklabel',yBinsStrT)
            set(gca,'xtick',[])
            set(gca,'fontsize',16)
            set(gca,'xtick',out.regionCutoffs);
            set(gca,'xticklabel','|');
            set(gca,'xticklabel',out.regions);
            set(gca,'XTickLabelRotation',270)
            xlabel('Electrodes','fontsize',20)
            ylabel('Time Bin','fontsize',20)
            print('-dpng','-loose',fname);
            
        end
        
        
        %% FIGURE 3 - time x freq (averaged across elecs)
        fname = fullfile(figDir,[subj 'Spect_TxF.png']);
        figs_subj.Spect_TxF = fname;
        spectTimeFreq(:,:,s)       = squeeze(mean(abs(out.meanWeightsPerTime),2));
        spectTimeFreqNonAbs(:,:,s) = squeeze(mean((out.meanWeightsPerTime),2));
        if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
            figure(3)
            clf
            imagesc(squeeze(mean(abs(out.meanWeightsPerTimeSort),2)));
            colormap('jet')
            colorbar
            
            set(gca,'xtick',1:nTimes)
            set(gca,'xticklabel',yBinsStrT)
            set(gca,'ytick',1:nFreqs)
            set(gca,'yticklabel',yBinsStrF)
            xlabel('Time Bin','fontsize',20)
            ylabel('Frequency','fontsize',20)
            set(gca,'fontsize',20)
            print('-dpng','-loose',fname);
        end
        
        %% FIGURE 4 - mean absolute weights by region for each time
        labels = params.timeBinLabels;
        if isempty(labels);labels=repmat({''},1,size(params.timeBins,1));end
        regions     = out.regionsAll;
%         regions = out.regionsAll;
%         regions = {'meanWeightsPerTimeHipp','meanWeightsPerTimeEC',...
%             'meanWeightsPerTimeMTL','meanWeightsPerTimeCA1',...
%             'meanWeightsPerTimeCA3','meanWeightsPerTimeDG',...
%             'meanWeightsPerTimeSub','meanWeightsPerTimePHC',...
%             'meanWeightsPerTimePRC','meanWeightsPerTimeTC',...
%             'meanWeightsPerTimeFC','meanWeightsPerTimeOC',...
%             'meanWeightsPerTimePC','meanWeightsPerTimeOth'};
        
        % this is stupid. I'm basically doing this twice so that I can find
        % out what the axis ranges will be and set them all the same
        ylims = [0 0];
        for t = 1:size(params.timeBins,1)
            regionData        = NaN(1,length(regions));
            regionDataNonZero = NaN(1,length(regions));
            for r = 1:length(regions)
                weights       = out.(regions{r})(:,:,t);
                isZero        = weights == 0;
                regionData(r)        = mean(nanmean(abs(weights),2));
                regionDataAll(r,s,t) = regionData(r);
                regionDataAllFreqsTimes(r,s,t,:) = nanmean(weights,2);                                
                
                regionDataNonZero(r)        = nanmean(abs(weights(~isZero)));
                regionDataNonZeroAll(r,s,t) = regionDataNonZero(r);
            end
            ylims(1) = max([ylims(1) max(regionData)]);
            ylims(2) = max([ylims(2) max(regionDataNonZero)]);
        end
        ylims = ceil(ylims * 100)/100;
        ylims(ylims==0) = .001;
        
        % loop over each time bin
        for t = 1:size(params.timeBins,1)
            
            % average weights within region, both including and excluding
            % zero weights
            regionData        = NaN(1,length(regions));
            regionDataNonZero = NaN(1,length(regions));
            for r = 1:length(regions)
                weights       = out.(regions{r})(:,:,t);
                isZero        = weights == 0;
                regionData(r)        = mean(nanmean(abs(weights),2));
                regionDataNonZero(r) = nanmean(abs(weights(~isZero)));
            end
            
            % plot bar
            plotData = {regionData,regionDataNonZero};
            ylabels  = {'Mean Abs Weights','Mean Abs NonZero Weights'};
            for panel = 1:2
                
                fname = fullfile(figDir,[subj '_bar_region_' labels{t} '_' strrep(ylabels{panel},' ','_') '.eps']);
                figs_subj.Region_Bar{t,panel} = fname;
                
                if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
                    figure(4)
                    clf
                    
                    %                     subplot(2,1,panel)
                    h=bar(plotData{panel},'linewidth',2,'facecolor',[.6 .6 .6]);
                    set(gca,'xtick',1:length(regions));
                    xNames = regions;
                    [xNames{isnan(plotData{panel})}] = deal('');
                    set(gca,'xticklabel',xNames)
                    ylabel(ylabels{panel});
                    set(gca,'fontsize',16)
                    set(gca,'xlim',[0 length(regions)+1])
                    set(gca,'ylim',[0 ylims(panel)])
                    grid on
                    set(gca,'gridlinestyle',':');
                    h=title([labels{t} ' Period'],'fontsize',20);
                    set(h,'fontweight','normal');
                    print('-depsc2','-loose',fname);
                end
                
            end
        end
        
        % also keep track of the weights from ther best time bin
        for r = 1:length(regions)
            weights = out.(regions{r})(:,:,out.bestTime);
            regionDataAllFreqsBestTime(r,s,:) = nanmean((weights),2);
        end
        
    end
    figs = [figs;figs_subj];
end

% also make group plots/report
fprintf('Creating group plots.\n');
figs_group = [];
figs_group.Region_Bar = {};


%% FIGURE - average time x freq plot
fname = fullfile(figDir,['group_Spect_TxF.png']);
figs_group.group_Spect_TxF = fname;
figs_group.N = sum(~isnan(spectTimeFreq(1,1,:)),3);
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(3)
    clf
    imagesc(nanmean(spectTimeFreq,3));
    colormap('jet')
    colorbar
    
    set(gca,'xtick',1:nTimes)
    set(gca,'xticklabel',yBinsStrT)
    set(gca,'ytick',1:nFreqs)
    set(gca,'yticklabel',yBinsStrF)
    xlabel('Time Bin','fontsize',20)
    ylabel('Frequency','fontsize',20)
    set(gca,'fontsize',20)
    print('-dpng','-loose',fname);
end


%% FIGURE - mean absolute weights for frequencies avg acros times
fname = fullfile(figDir,['group_freq_bar.eps']);
figs_group.group_freq_bar = fname;

% figs_group.N = sum(~isnan(spectTimeFreq(1,1,:)),3);
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(3)
    clf
    
    plotData = squeeze(nanmean(spectTimeFreq(:,1:7,:),2));
    err      = nanstd(plotData,[],2)/sqrt(figs_group.N-1);
    h=bar(nanmean(plotData,2),'linewidth',2,'facecolor',[.6 .6 .6]);
    hold on
    errorbar(1:size(plotData,1),nanmean(plotData,2),err,'k','linewidth',2,'linestyle','none')
    
    set(gca,'xtick',1:nFreqs)
    set(gca,'xticklabel',yBinsStrF)
    ylabel('Mean Abs Weights','fontsize',20)
    xlabel('Frequency','fontsize',20)
    set(gca,'fontsize',20)
    grid on
    set(gca,'xlim',[0 nFreqs+1])
    set(gca,'gridlinestyle',':');
    print('-depsc2','-loose',fname);
end

%% FIGURE - mean absolute weights for times avg acros freqs
fname = fullfile(figDir,['group_time_bar.eps']);
figs_group.group_time_bar = fname;
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(3)
    clf
    
    plotData = squeeze(nanmean(spectTimeFreq,1));
    err      = nanstd(plotData,[],2)/sqrt(figs_group.N-1);
    h=bar(nanmean(plotData,2),'linewidth',2,'facecolor',[.6 .6 .6]);
    hold on
    errorbar(1:size(plotData,1),nanmean(plotData,2),err,'k','linewidth',2,'linestyle','none')
    
    set(gca,'xtick',1:nTimes)
    set(gca,'xticklabel',yBinsStrT)
    ylabel('Mean Abs Weights','fontsize',20)
    xlabel('Time Bin','fontsize',20)
    set(gca,'fontsize',20)
    grid on
    set(gca,'gridlinestyle',':');
    set(gca,'xlim',[0 nTimes+1])
    print('-depsc2','-loose',fname);
end

%% FIGURE - mean absolute weights for times avg acros freqs
fname = fullfile(figDir,['group_region_bar.eps']);
figs_group.group_region_bar = fname;
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(3)
    clf
    
    plotData = nanmean(nanmean(regionDataAll,3),2);
    err      = nanstd(nanmean(regionDataAll,3),[],2)./sqrt(sum(~isnan(nanmean(regionDataAll,3)),2)-1);
    h=bar(nanmean(plotData,2),'linewidth',2,'facecolor',[.6 .6 .6]);
    hold on
    errorbar(1:size(plotData,1),nanmean(plotData,2),err,'k','linewidth',2,'linestyle','none')
    set(gca,'xtick',1:length(regions));
    xNames = regions;
    %             [xNames{isnan(plotData{panel})}] = deal('');
    set(gca,'xticklabel',xNames)
    ylabel(ylabels{1});
    set(gca,'fontsize',16)
    set(gca,'xlim',[0 length(regions)+1])
    %             set(gca,'ylim',[0 ylims(panel)])
    grid on
    set(gca,'gridlinestyle',':');
    
    print('-depsc2','-loose',fname);
end

%% FIGURE - mean absolute weights for times avg acros freqs
fname = fullfile(figDir,['group_regionNonZero_bar.eps']);
figs_group.group_regionNonZero_bar = fname;
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(3)
    clf
    
    plotData = nanmean(nanmean(regionDataNonZeroAll,3),2);
    err      = nanstd(nanmean(regionDataNonZeroAll,3),[],2)./sqrt(sum(~isnan(nanmean(regionDataNonZeroAll,3)),2)-1);
    h=bar(nanmean(plotData,2),'linewidth',2,'facecolor',[.6 .6 .6]);
    hold on
    errorbar(1:size(plotData,1),nanmean(plotData,2),err,'k','linewidth',2,'linestyle','none')
    set(gca,'xtick',1:length(regions));
    xNames = regions;
    %             [xNames{isnan(plotData{panel})}] = deal('');
    set(gca,'xticklabel',xNames)
    ylabel(ylabels{2});
    set(gca,'fontsize',16)
    set(gca,'xlim',[0 length(regions)+1])
    %             set(gca,'ylim',[0 ylims(panel)])
    grid on
    set(gca,'gridlinestyle',':');
    print('-depsc2','-loose',fname);
end

%% FIGURE - mean absolute weights by region for each time
for t = 1:size(params.timeBins,1)
    
    % plot bar
    clf
    err1 =  nanstd(regionDataAll(:,:,t),[],2)./sqrt(sum(~isnan(regionDataAll(:,:,t)),2)-1);
    err2 =  nanstd(regionDataNonZeroAll(:,:,t),[],2)./sqrt(sum(~isnan(regionDataNonZeroAll(:,:,t)),2)-1);
    plotDataY = {nanmean(regionDataAll(:,:,t),2),nanmean(regionDataNonZeroAll(:,:,t),2)};
    errDataY  = {err1 err2};
    ylabels  = {'Mean Abs Weights','Mean Abs NonZero Weights'};
    for panel = 1:2
        
        fname = fullfile(figDir,['group_bar_region_' labels{t} '_' strrep(ylabels{panel},' ','_') '.eps']);
        figs_group.Region_Bar{t,panel} = fname;
        
        if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
            figure(4)
            clf
            h=bar(plotDataY{panel},'linewidth',2,'facecolor',[.6 .6 .6]);
            hold on
            errorbar(1:length(regions),plotDataY{panel},errDataY{panel},'k','linewidth',2,'linestyle','none')
            set(gca,'xtick',1:length(regions));
            xNames = regions;
            %             [xNames{isnan(plotData{panel})}] = deal('');
            set(gca,'xticklabel',xNames)
            ylabel(ylabels{panel});
            set(gca,'fontsize',16)
            set(gca,'xlim',[0 length(regions)+1])
            %             set(gca,'ylim',[0 ylims(panel)])
            grid on
            set(gca,'gridlinestyle',':');
            h=title([labels{t} ' Period'],'fontsize',20);
            set(h,'fontweight','normal');
            print('-depsc2','-loose',fname);
        end
    end
end
keyboard
%% FIGURE - mean (non-abs) weights by region for best time for each freq
mWeights = squeeze(nanmean(regionDataAllFreqsBestTime,2));
sWeights = squeeze(nanstd(regionDataAllFreqsBestTime,[],2));
eWeights = sWeights./sqrt(squeeze(sum(~isnan(regionDataAllFreqsBestTime),2))-1);
[h,p,c,s] = ttest(permute(regionDataAllFreqsBestTime,[2 1 3]));
h=squeeze(h);p=squeeze(p);t=squeeze(s.tstat);

for f = 1:size(mWeights,2)

    figure(5)
    clf
    hold on
    for r = 1:size(mWeights,1)
        c = [.7 .7 .7];
        if ~isnan(h(r,f)) && h(r,f)
            c = [.2 .2 .7];
            if t(r,f) > 0
                c = [.7 .2 .2];
            end
        end                
        bar(r,mWeights(r,f),'linewidth',2,'FaceColor',c)
    end
    errorbar(1:r,mWeights(:,f),eWeights(:,f),'k','linewidth',2,'linestyle','none')
    set(gca,'xtick',1:length(regions));
    set(gca,'xticklabel',regions)
    set(gca,'XTickLabelRotation',45)
    xlabel('Region','fontsize',16)
    ylabel('Average Classifier Weight','fontsize',16)
    set(gca,'fontsize',16)
    set(gca,'xlim',[0 length(regions)+1]);
    grid on
    set(gca,'gridlinestyle',':')
    
    fname = fullfile(figDir,['group_bar_region_freq_' num2str(f) '.eps']);
    print('-depsc2','-loose',fname); 
end

%% FIGURE - mean (non-abs) weights by region for best time for binned freqs
%%% HARD CODED STUFF FIX IT
lf = 1:4;
hf = 8:11;
regionDataLFHF = NaN(9,64,2);
regionDataLFHF(:,:,1) = nanmean(regionDataAllFreqsBestTime(1:9,:,lf),3);
regionDataLFHF(:,:,2) = nanmean(regionDataAllFreqsBestTime(1:9,:,hf),3);
mWeights = squeeze(nanmean(regionDataLFHF,2));            
sWeights = squeeze(nanstd(regionDataLFHF,[],2));            
eWeights = sWeights./sqrt(squeeze(sum(~isnan(regionDataLFHF),2))-1);            
[h,p,c,s] = ttest(permute(regionDataLFHF,[2 1 3]));            
h=squeeze(h);p=squeeze(p);t=squeeze(s.tstat);



for f = 1:size(mWeights,2)

    figure(5)
    clf
    hold on
    for r = 1:size(mWeights,1)
        c = [.7 .7 .7];
        if ~isnan(h(r,f)) && h(r,f)
            c = [.2 .2 .7];
            if t(r,f) > 0
                c = [.7 .2 .2];
            end
        end                
        bar(r,mWeights(r,f),'linewidth',2,'FaceColor',c)
    end
    errorbar(1:r,mWeights(:,f),eWeights(:,f),'k','linewidth',2,'linestyle','none')
    set(gca,'xtick',1:length(regions)-1);
%     set(gca,'xticklabel',regions(1:end-1))
    set(gca,'XTickLabelRotation',45)
%     xlabel('Region','fontsize',16)
    ylabel('Average Classifier Weight','fontsize',16)
    set(gca,'fontsize',16)
    set(gca,'xlim',[0 length(regions)]);
    grid on
    set(gca,'gridlinestyle',':')
    
    fname = fullfile(figDir,['group_bar_region_freqBin_' num2str(f) '.eps']);
    print('-depsc2','-loose',fname); 
end

%% FIGURE - mean (non-abs) weights time by freq spectrogram
figure(6)
clf

% make x labels
xBins    = round((params.timeBins) / 10) / 100;
f = @(x) [num2str(x(1)), '-', num2str(x(2))];
xlabels = cellfun(f,num2cell(xBins,2),'uniformoutput',false);

%%%% GET THIS FROM PARAMS
% make y labels
freqs = logspace(log10(3),log10(120),12);
f = @(x) sprintf('%.1f',x);
ylabels = cellfun(f,num2cell(freqs),'uniformoutput',false);

% plot it
mTF = flipud(nanmean(spectTimeFreqNonAbs,3));
imagesc(mTF);
set(gca,'xtick',1:10:size(mTF,2));
set(gca,'xticklabel',xlabels(1:10:size(mTF,2)))
clim = get(gca,'clim');
set(gca,'clim',[-abs(max(clim)) abs(max(clim))])
colorbar

% set labels
xlabel('Time (s)','fontsize',20)
ylabel('Frequency','fontsize',20)
set(gca,'ytick',1:length(freqs))
set(gca,'yticklabel',fliplr(ylabels))
set(gca,'fontsize',20)
set(gca,'xticklabelrotation',45)
colormap default

fname = fullfile(figDir,['TF_nonAbs.png']);
print('-dpng','-loose',fname);

%%%%%%%%%%%%%%%%%



good = ~cellfun('isempty',{figs.subj});
figs = figs(good);
texName = 'lassoWeights_report.tex';
write_texfile(saveDir,texName,figs)


curr_dir = pwd;
cd(saveDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(saveDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);

% compile group report
texName = 'group_lassoWeights_report.tex';
write_texfile_group(saveDir,texName,figs_group)


curr_dir = pwd;
cd(saveDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(saveDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);
keyboard


function out = loadWeightsByRegion(subj,saveDir)
out = [];

lassoFile  = fullfile(saveDir,[subj '_lasso.mat']);
if ~exist(lassoFile,'file')
    fprintf('No lasso for %s.\n',subj)
    out = [];
    return
end

% load classification model
lassoData  = load(lassoFile);

% get subject electrode info
tal                   = lassoData.tal;

% number of features
nFeatures = length(lassoData.res(1).A{1});
nFreqs    = size(lassoData.params.freqBins,1);
nTimes    = size(lassoData.params.timeBins,1);

% bestTime      = find(lassoData.AUC == max(lassoData.AUC),1,'first');
bestTime = 9;
out.bestTime  = bestTime;

if lassoData.params.modelEachTime
    nElecs = nFeatures/nFreqs;
else
    nElecs = nFeatures/nFreqs/nTimes;
end

if nElecs ~= length(tal)
    fprintf('Number of electrodes in tal structure does not match features.\n This should not be possible.\n')
    return
end

% make sure loctag isn't missing
if ~isfield(tal,'locTag')
    [tal.locTag] = deal('');
end
if sum(cellfun('isempty',{tal.locTag})) == length(tal)
    [tal.locTag] = deal('');
end
missing               = cellfun('isempty',{tal.locTag});
[tal(missing).locTag] = deal('');


% get the electrode indices of brain regions
% locTag based
elecs = [];
elecs.H       = ~cellfun('isempty',regexpi({tal.locTag},['CA1|CA2|CA3|DG|sub']));
elecs.ec      = ~cellfun('isempty',regexpi({tal.locTag},['ec|erc']));
elecs.MTL     = ~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA2|CA3|DG|sub|amy|phc|prc|BA36|erc']));
elecs.MTL     = ~cellfun('isempty',regexpi({tal.locTag},['ec|amy|phc|prc|BA36|erc']));
elecs.ca1     = ~cellfun('isempty',regexpi({tal.locTag},['ca1']));
elecs.ca3     = ~cellfun('isempty',regexpi({tal.locTag},['ca3']));
elecs.dg      = ~cellfun('isempty',regexpi({tal.locTag},['dg']));
elecs.sub     = ~cellfun('isempty',regexpi({tal.locTag},['sub']));
elecs.phc     = ~cellfun('isempty',regexpi({tal.locTag},['phc']));
elecs.prc     = ~cellfun('isempty',regexpi({tal.locTag},['prc']));

% lobe based
elecs.frontal = strcmp({tal.Loc2},'Frontal Lobe');
elecs.occ     = strcmp({tal.Loc2},'Occipital Lobe');
elecs.par     = strcmp({tal.Loc2},'Parietal Lobe');
elecs.temp    = strcmp({tal.Loc2},'Temporal Lobe') & ~elecs.MTL & ~elecs.H;

% new version based on brodmann areas
ba = {tal.Loc5};
%%%% THE SPACE BEFORE THE NUMBER IS IMPORTANT %%%
elecs.aPFC = ~cellfun('isempty',regexpi(ba,[' 10| 11'])) & ~elecs.MTL & ~elecs.H;
elecs.mPFC = ~cellfun('isempty',regexpi(ba,[' 24| 25| 32| 33'])) & ~elecs.MTL & ~elecs.H;
elecs.PFC  = ~cellfun('isempty',regexpi(ba,[' 45| 47| 9| 46'])) & ~elecs.MTL & ~elecs.H;
elecs.TC   = ~cellfun('isempty',regexpi(ba,[' 20| 21| 37'])) & ~elecs.MTL & ~elecs.H;
elecs.PPC  = ~cellfun('isempty',regexpi(ba,[' 7| 40| 39'])) & ~elecs.MTL & ~elecs.H;
elecs.mPC  = ~cellfun('isempty',regexpi(ba,[' 23| 29| 30| 31'])) & ~elecs.MTL & ~elecs.H;
elecs.OC   = ~cellfun('isempty',regexpi(ba,[' 17| 18| 19'])) & ~elecs.MTL & ~elecs.H;

% average across all folds and reshape into freq x elec x time.
meanWeightsPerTime = NaN(nFreqs,nElecs,nTimes);
if lassoData.params.modelEachTime
    for t = 1:length(lassoData.res)
        for fold = 1:length(lassoData.res(t).A)
            if isrow(lassoData.res(t).A{fold})
                lassoData.res(t).A{fold} = lassoData.res(t).A{fold}';
            end
        end
        meanTmp = mean(horzcat(lassoData.res(t).A{:}),2);
        meanTmp = reshape(meanTmp,nFreqs,nElecs);
        meanWeightsPerTime(:,:,t) = meanTmp;
    end
else
    for fold = 1:length(lassoData.res.A)
        if isrow(lassoData.res.A{fold})
            lassoData.res.A{fold} = lassoData.res.A{fold}';
        end
    end
    meanTmp = mean(horzcat(lassoData.res.A{:}),2);
    meanWeightsPerTime = reshape(meanTmp,nFreqs,nElecs,[]);
end

% filter by regions
regions = {'aPFC','mPFC','PFC','MTL','H','TC','PPC','mPC','OC'};

% all elecs
out.meanWeightsPerTime     = meanWeightsPerTime;

% elecs by region
elecs.Other = true(1,length(tal));
elec_order = [];
regionCount = zeros(1,length(regions)+1);
for r = 1:length(regions)
    out.(regions{r}) = meanWeightsPerTime(:,elecs.(regions{r}),:);  
    elec_order = [elec_order find(elecs.(regions{r}))];
    elecs.Other = elecs.Other & ~elecs.(regions{r});
    regionCount(r) = sum(elecs.(regions{r}));
end
out.Other = meanWeightsPerTime(:,elecs.Other,:); 
elec_order = [elec_order find(elecs.Other)];
regions{end+1} = 'Other';
regionCount(end) = sum(elecs.Other);

% created a sorted by region version
out.meanWeightsPerTimeSort = meanWeightsPerTime(:,elec_order,:);
noElecs = regionCount == 0;
out.regionCutoffs = cumsum(regionCount(~noElecs));
out.regions       = regions(~noElecs);
out.regionsAll    = regions;



% Start making the tex file
function write_texfile(saveDir,texName, figs)

% Write the document. If you do not have write permission, this will crash.
fid = fopen(fullfile(saveDir,texName),'w');

if fid==-1;
    error(sprintf('cannot open %s',texName))
end

% Write out the preamble to the tex doc. This is standard stuff and doesn't
% need to be changed
fprintf(fid,'\\documentclass[a4paper]{article} \n');
fprintf(fid,'\\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}\n');
fprintf(fid,'\\usepackage{graphicx,multirow} \n');
fprintf(fid,'\\usepackage{epstopdf} \n');
fprintf(fid,'\\usepackage[small,bf,it]{caption}\n');
fprintf(fid,'\\usepackage{subfig,amsmath} \n');
fprintf(fid,'\\usepackage{wrapfig} \n');
fprintf(fid,'\\usepackage{longtable} \n');
fprintf(fid,'\\usepackage{pdfpages}\n');
fprintf(fid,'\\usepackage{mathtools}\n');
fprintf(fid,'\\usepackage{array}\n');
fprintf(fid,'\\usepackage{enumitem}\n');
fprintf(fid,'\\usepackage{sidecap} \\usepackage{soul}\n');

% fprintf(fid,'\\setlength\\belowcaptionskip{5pt}\n');
fprintf(fid,'\n');
fprintf(fid,'\\addtolength{\\oddsidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\evensidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\textwidth}{1.75in} \n');
fprintf(fid,'\\addtolength{\\topmargin}{-.75in} \n');
fprintf(fid,'\\addtolength{\\textheight}{1.75in} \n');
fprintf(fid,'\n');
fprintf(fid,'\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}} \n');

fprintf(fid,'\\usepackage{fancyhdr}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
% fprintf(fid,'\\lhead{Report: %s }\n',strrep(subj,'_','\_'));
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
for s = 1:length(figs)
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\includegraphics[width=0.6\\textwidth]{%s}\n',figs(s).Spect_Freq);
    fprintf(fid,'\\includegraphics[width=0.6\\textwidth]{%s}\n',figs(s).Spect_Time);
    fprintf(fid,'\\includegraphics[width=0.6\\textwidth]{%s}\n',figs(s).Spect_TxF);
    fprintf(fid,'\\caption{%s. Absolute classifier weights (including zeros).}\n\n',strrep(figs(s).subj,'_',' '));
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    for r = 1:size(figs(s).Region_Bar,1)
        fprintf(fid,'\\includegraphics[width=0.3\\textwidth]{%s}\n',figs(s).Region_Bar{r,1});
    end
    fprintf(fid,'\\caption{Average weights across electrodes by brain region, including zeros, for each time bin.}\n\n');
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    for r = 1:size(figs(s).Region_Bar,1)
        fprintf(fid,'\\includegraphics[width=0.3\\textwidth]{%s}\n',figs(s).Region_Bar{r,2});
    end
    fprintf(fid,'\\caption{Average weights across electrodes by brain region, excluding zeros, for each time bin.}\n\n');
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\clearpage\n\n\n');
    
end

fprintf(fid,'\\end{document}\n\n\n');












% Start making the tex file
function write_texfile_group(saveDir,texName, figs)

% Write the document. If you do not have write permission, this will crash.
fid = fopen(fullfile(saveDir,texName),'w');

if fid==-1;
    error(sprintf('cannot open %s',texName))
end

% Write out the preamble to the tex doc. This is standard stuff and doesn't
% need to be changed
fprintf(fid,'\\documentclass[a4paper]{article} \n');
fprintf(fid,'\\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}\n');
fprintf(fid,'\\usepackage{graphicx,multirow} \n');
fprintf(fid,'\\usepackage{epstopdf} \n');
fprintf(fid,'\\usepackage[small,bf,it]{caption}\n');
fprintf(fid,'\\usepackage{subfig,amsmath} \n');
fprintf(fid,'\\usepackage{wrapfig} \n');
fprintf(fid,'\\usepackage{longtable} \n');
fprintf(fid,'\\usepackage{pdfpages}\n');
fprintf(fid,'\\usepackage{mathtools}\n');
fprintf(fid,'\\usepackage{array}\n');
fprintf(fid,'\\usepackage{enumitem}\n');
fprintf(fid,'\\usepackage{sidecap} \\usepackage{soul}\n');

% fprintf(fid,'\\setlength\\belowcaptionskip{5pt}\n');
fprintf(fid,'\n');
fprintf(fid,'\\addtolength{\\oddsidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\evensidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\textwidth}{1.75in} \n');
fprintf(fid,'\\addtolength{\\topmargin}{-.75in} \n');
fprintf(fid,'\\addtolength{\\textheight}{1.75in} \n');
fprintf(fid,'\n');
fprintf(fid,'\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}} \n');

fprintf(fid,'\\usepackage{fancyhdr}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
% fprintf(fid,'\\lhead{Report: %s }\n',strrep(subj,'_','\_'));
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
for s = 1:length(figs)
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\includegraphics[width=0.6\\textwidth]{%s}\n',figs(s).group_Spect_TxF);
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).group_freq_bar);
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).group_time_bar);
    fprintf(fid,'\\caption{%d subjects. Top: Time x Frequency spectrogram of subject average absolute classifier weights (includes zeros). Bottom Left: Averaged across time. Bottom Right: Averaged across frequncies.}\n\n',figs.N);
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).group_region_bar);
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).group_regionNonZero_bar);
    fprintf(fid,'\\caption{%d subjects. Average weights across subjects by brain region. Left: Including zero weights. Right: Exlcuding zero weights.}\n\n',figs.N);
    fprintf(fid,'\\end{figure}\n\n\n');
    fprintf(fid,'\\clearpage\n\n\n');
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    for r = 1:size(figs(s).Region_Bar,1)
        fprintf(fid,'\\includegraphics[width=0.28\\textwidth]{%s}\n',figs(s).Region_Bar{r,1});
    end
    fprintf(fid,'\\caption{Average weights across subjects by brain region, including zeros, for each time bin.}\n\n');
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    for r = 1:size(figs(s).Region_Bar,1)
        fprintf(fid,'\\includegraphics[width=0.28\\textwidth]{%s}\n',figs(s).Region_Bar{r,2});
    end
    fprintf(fid,'\\caption{Average weights across subjects by brain region, excluding zeros, for each time bin.}\n\n');
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\clearpage\n\n\n');
    
end

fprintf(fid,'\\end{document}\n\n\n');






