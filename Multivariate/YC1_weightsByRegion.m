function YC1_weightsByRegion(subjs,params)

if ~exist('params','var') || isempty(params)
    params = multiParams();
end

overwrite = true

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

subjs = subjs(~strcmp(subjs,'R1006P'));
subjs = subjs(~strcmp(subjs,'R1059J_1'));
subjs = subjs(~strcmp(subjs,'R1061T'));
keyboard

nTimes = size(params.timeBins,1);
nFreqs = size(params.freqBins,1);
spectTimeFreq        = NaN(size(params.freqBins,1),size(params.timeBins,1),length(subjs));
regionDataAll        = NaN(8,length(subjs),size(params.timeBins,1));
regionDataNonZeroAll = NaN(8,length(subjs),size(params.timeBins,1));
movementClass        = {};

figs = [];
for s = 1:length(subjs)    
    subj = subjs{s};
    fprintf('Processing %s.\n',subjs{s})
    
    figs_subj = struct('subj',[],'Spect_Freq',[],'Spect_Time',[],...
        'Spect_TxF',[],'Region_Bar',[],'nElecs',[]);
    
    % load the weights for this subject
    out = loadWeightsByRegion(subjs{s},dataDir);
    if ~isempty(out)        
        figs_subj.subj = subj;
        figs_subj.Region_Bar = {};
        movementClass{s} = out.movementClass;
        
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
        spectTimeFreq(:,:,s) = squeeze(mean(abs(out.meanWeightsPerTimeSort),2));
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
        labels = {'Pre','Spin','Drive','Wait','Post'};
        if size(params.timeBins,1) == 1
            labels = {'Encoding'};
        end
        if size(params.timeBins,1) == 6
            labels = {'Pre','Spin','Drive','Wait','Post','Enc'};
        end
        if size(params.timeBins,1) == 8
            labels = {'Pre','Spin','Drive','Drive','Drive','Wait','Post','Enc'};
        end        
        regions     = {'H','EC','MTL','TC','FC','OC','PC','X'};
        fieldsToUse = {'meanWeightsPerTimeHipp','meanWeightsPerTimeEC',...
            'meanWeightsPerTimeMTL','meanWeightsPerTimeTC',...
            'meanWeightsPerTimeFC','meanWeightsPerTimeOC',...
            'meanWeightsPerTimePC','meanWeightsPerTimeOth'};
        
        % this is stupid         
        ylims = [0 0];
        for t = 1:size(params.timeBins,1)
            regionData        = NaN(1,length(fieldsToUse));
            regionDataNonZero = NaN(1,length(fieldsToUse));            
            for r = 1:length(fieldsToUse)
                weights       = out.(fieldsToUse{r})(:,:,t);
                isZero        = weights == 0;
                regionData(r)        = mean(nanmean(abs(weights),2));
                regionDataAll(r,s,t) = regionData(r);
                
                regionDataNonZero(r)        = nanmean(abs(weights(~isZero)));
                regionDataNonZeroAll(r,s,t) = regionDataNonZero(r);
            end            
            ylims(1) = max([ylims(1) max(regionData)]);
            ylims(2) = max([ylims(2) max(regionDataNonZero)]);            
        end
        ylims = ceil(ylims * 100)/100;    
        
        % loop over each time bin
        for t = 1:size(params.timeBins,1)                        
            
            % average weights within region, both including and excluding
            % zero weights
            regionData        = NaN(1,length(fieldsToUse));
            regionDataNonZero = NaN(1,length(fieldsToUse));
            for r = 1:length(fieldsToUse)
                weights       = out.(fieldsToUse{r})(:,:,t);
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
    
    plotData = squeeze(nanmean(spectTimeFreq,2));
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
% if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
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
%     h=title([labels{t} ' Period'],'fontsize',20);
%     set(h,'fontweight','normal');
    print('-depsc2','-loose',fname);    
% end

%% FIGURE - mean absolute weights for times avg acros freqs
fname = fullfile(figDir,['group_regionNonZero_bar.eps']);
figs_group.group_regionNonZero_bar = fname;
% if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
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
%     h=title([labels{t} ' Period'],'fontsize',20);
%     set(h,'fontweight','normal');
    print('-depsc2','-loose',fname);    
% end

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
%%%%%%%%%%%%%%%%%

keyboard

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

chanceFile = fullfile(saveDir,[subj '_chance_perf_dist.mat']);
lassoFile  = fullfile(saveDir,[subj '_lasso.mat']);
if ~exist(lassoFile,'file') || ~exist(chanceFile,'file')
    fprintf('No lasso/chance file for %s.\n',subj)
    out = [];
    return
end

% load classification model
lassoData  = load(lassoFile);
chanceData = load(chanceFile);

% get the significant (percetile) at each time point
perc = mean(repmat(lassoData.perf,size(chanceData.perf_all,1),1) > chanceData.perf_all);
[~,bestTimeBin] = max(perc);

% get subject electrode info
tal = lassoData.tal;

nFeatures = length(lassoData.res(1).A{1});
nFreqs    = size(lassoData.params.freqBins,1);
nTimes    = size(lassoData.params.timeBins,1);
if lassoData.params.modelEachTime
    nElecs = nFeatures/nFreqs;
else
    nElecs = nFeatures/nFreqs/nTimes;
end

if nElecs ~= length(tal)
    fprintf('Number of electrodes in tal structure does not match features.\n This should not be possible.\n')
    return
end

if ~isfield(tal,'locTag')
    fprintf('No loc tag information for %s.\n',subj)
    return
end

% get the electrode indices of brain regions
hipp_elecs    = ~cellfun('isempty',regexpi({tal.locTag},['CA1|CA3|DG|sub']));
ec_elecs      = ~cellfun('isempty',regexpi({tal.locTag},['ec|erc']));
mtl_elecs     = ~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc']));
frontal_elecs = strcmp({tal.Loc2},'Frontal Lobe');
occ_elecs     = strcmp({tal.Loc2},'Occipital Lobe');
par_elecs     = strcmp({tal.Loc2},'Parietal Lobe');
temp_elecs    = strcmp({tal.Loc2},'Temporal Lobe') & ~mtl_elecs;
other_elecs   = ~(hipp_elecs | ec_elecs | mtl_elecs | frontal_elecs | occ_elecs | par_elecs | temp_elecs);

% pval of different time bins
p = 1-perc;
thresh = p < .05;
sigTimes = p < thresh;
out.sigTimes = sigTimes;

% make this part of the params in the future
params.timeIsMovement = [false true true true true false false false];

if any(sigTimes(1:end-1)) && ~any(params.timeIsMovement(sigTimes(1:end-1))==0) && sum(sigTimes(1:end-1))>=2
    out.movementClass = 'mover';
elseif any(sigTimes(1:end-1)) && ~any(params.timeIsMovement(sigTimes(1:end-1))) && sum(sigTimes(1:end-1))>=2
    out.movementClass = 'nonMover';
elseif any(sigTimes(1:end-1)) && sum(sigTimes(1:end-1)) > 3
    out.movementClass = 'both';
else
    out.movementClass = 'neither';
end


out.meanAccNonmove = mean(lassoData.perf(params.timeIsMovement));
out.meanAccMove = mean(lassoData.perf(~params.timeIsMovement));

% average across all folds and reshape into freq x elec x time
meanWeightsPerTime = NaN(nFreqs,nElecs,nTimes);
for t = 1:length(lassoData.res)
    meanTmp = mean(horzcat(lassoData.res(t).A{:}),2);
    meanTmp = reshape(meanTmp,nFreqs,nElecs);
    meanWeightsPerTime(:,:,t) = meanTmp;
end

% filter by regions
out.meanWeightsPerTime     = meanWeightsPerTime;
out.meanWeightsPerTimeHipp = meanWeightsPerTime(:,hipp_elecs,:);
out.meanWeightsPerTimeEC   = meanWeightsPerTime(:,ec_elecs,:);
out.meanWeightsPerTimeMTL  = meanWeightsPerTime(:,mtl_elecs,:);
out.meanWeightsPerTimeFC   = meanWeightsPerTime(:,frontal_elecs,:);
out.meanWeightsPerTimeOC   = meanWeightsPerTime(:,occ_elecs,:);
out.meanWeightsPerTimePC   = meanWeightsPerTime(:,par_elecs,:);
out.meanWeightsPerTimeTC   = meanWeightsPerTime(:,temp_elecs,:);
out.meanWeightsPerTimeOth  = meanWeightsPerTime(:,other_elecs,:);

% [
elec_order = [find(hipp_elecs) find(ec_elecs) find(mtl_elecs) ...
              find(temp_elecs) find(frontal_elecs) find(occ_elecs) ...
              find(par_elecs) find(other_elecs)];
out.meanWeightsPerTimeSort = meanWeightsPerTime(:,elec_order,:);

regions     = {'H','EC','MTL','TC','FC','OC','PC','X'};
regionCount = sum([hipp_elecs' ec_elecs' ...
    mtl_elecs' temp_elecs' frontal_elecs'...
    occ_elecs' par_elecs' other_elecs']);
noElecs     = regionCount == 0;


out.regionCutoffs = cumsum(regionCount(~noElecs));
out.regions       = regions(~noElecs);


% out.regionCutoffs = cumsum(sum([hipp_elecs' ec_elecs' ...
%     mtl_elecs' temp_elecs' frontal_elecs'...
%     occ_elecs' par_elecs' other_elecs']));



           
% for r = 1:length(regions)
%     fprintf('%s: %.4f, %.4f percent\n',regions{r},out.wRegions(r),out.pRegions(r)*100)
% end




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






