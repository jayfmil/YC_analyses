function YC2_makeSubjectReports(subjs,params,overwrite)
% function YC2_makeSubjectReports(subjs,params,overwrite)
%
% Inputs:
%
%      subjs - cell array of subject strings (default get_subs('RAM_YC1'))
%     params - params structure (default is returned by multiParams)
%  overwrite - boolean. if true, overwrite existing figures
%
% Make a report of the classifier performance for each subject in YC1 that
% has classification and chance clasification data. Also make a group
% average report.
%
% For each subject, figures are:
%
%   - Classifier performance over time for % correct and for AUC
%   - Recall change as a function of classifier output for the best
%     performing  time bin
%   - Histogram of the patient's behavioral perfomance
%
%
% For the group report, figures are quartile plots for each time bin,
% accuracy historograms on the subject level for each time bin, and AUC
% histograms on the subject level for each time bin.
%
%  Reports are saved to <location of subject
%  data>/reports/lassoChance_report.pdf and group_lassoChance_report.pdf

% if not given, use default params
if ~exist('params','var') || isempty(params)
    params = multiParams();
end

% overwrite exisiting figures?
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

% tex directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
dataDir = fullfile(params.basePath,f(params.useCorrectedPower,y),'YC2');
saveDir = fullfile(dataDir,'reports');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC2');
end
% subjs = subjs(~strcmp(subjs,'R1061T'))

% figure directory
figDir = fullfile(saveDir,'figs');
if ~exist('figDir','dir')
    mkdir(figDir)
end

% store group information
perf_all      = NaN(length(subjs),size(params.timeBins,1));
perf_p_all    = NaN(length(subjs),size(params.timeBins,1));
auc_all       = NaN(length(subjs),size(params.timeBins,1));
aucNull_all   = NaN(length(subjs),size(params.timeBins,1));
aucNull_all   = [];
auc_p_all     = NaN(length(subjs),size(params.timeBins,1));
perfEnc_all   = NaN(length(subjs),size(params.timeBins,1));
perfEnc_p_all = NaN(length(subjs),size(params.timeBins,1));
aucEnc_all    = NaN(length(subjs),size(params.timeBins,1));
aucEnc_p_all  = NaN(length(subjs),size(params.timeBins,1));
quarts_all    = NaN(size(params.timeBins,1),3,length(subjs));
quartsEnc_all = NaN(size(params.timeBins,1),3,length(subjs));
quartsBest_all = NaN(size(params.timeBins,1),3,length(subjs));
quartsBestEnc_all = NaN(size(params.timeBins,1),3,length(subjs));

aucBest_all       = NaN(length(subjs),size(params.timeBins,1));
aucBest_p_all     = NaN(length(subjs),size(params.timeBins,1));

aucBestEnc_all    = NaN(length(subjs),size(params.timeBins,1));
aucBestEnc_p_all  = NaN(length(subjs),size(params.timeBins,1));
% will hold figure paths for latex report
figs = [];
for s = 1:length(subjs)
    subj = subjs{s};
    fprintf('Creating plots for %s.\n',subj);
    
    % will subject specific figure paths
    figs_subj = struct('subj',[],'behavior',[],'quartsperf',[],'region','',...
        'AUC','','perf','','nElecs',[],'quartsAUC',[]);
    if isempty(params.region)
        params.region = 'all';
    end
    figs_subj.region = params.region;
    figs_subj.subj   = subj;
    
    % see if files exist for subject. if not, continue
    chanceFile = fullfile(dataDir,[subj '_YC2_chance_perf_dist.mat']);
    lassoFile  = fullfile(dataDir,[subj '_YC2_lasso.mat']);
    if ~exist(chanceFile,'file')
        fprintf('Chance distribution not found for %s.\n',subj)
        continue
    end
    if ~exist(lassoFile,'file')
        fprintf('Lasso not found for %s.\n',subj)
        continue
    end
    
    events     = get_sub_events('RAM_YC2',subj);
    
    % load lasso and chance lasso data
    chanceData = load(chanceFile);
    lassoData  = load(lassoFile);
    figs_subj.nElecs = length(lassoData.tal);
    
    % 
%     p = mean(repmat(lassoData.AUC,[size(chanceData.auc_all,1), 1]) > chanceData.auc_all);
%     maxP = max(p);
%     AUC = yc1Data.AUC;
%     AUC(p~=maxP) = NaN;
%     [~,timeToUse] = max(AUC);    
%     keyboard
    % calculate percentile for subject for the timebin to timebin direct
    % mapping
    perf_p              = mean(repmat(lassoData.perf,size(chanceData.perf_all,1),1) > chanceData.perf_all);
    auc_p               = mean(repmat(lassoData.AUC,size(chanceData.auc_all,1),1) > chanceData.auc_all);
    aucBest_p           = mean(repmat(lassoData.AUC_best,size(chanceData.aucBest_all,1),1) > chanceData.aucBest_all);
    perf_p_all(s,:)     = perf_p;
    auc_p_all(s,:)      = auc_p;
    aucBest_p_all(s,:)  = aucBest_p;
    perf_all(s,:)       = lassoData.perf;
    auc_all(s,:)        = lassoData.AUC;
    aucNull_all    = [aucNull_all; chanceData.auc_all];
    aucBest_all(s,:)    = lassoData.AUC_best;
    
    % now for timebin to average encoding interval
    perf_pEnc              = mean(repmat([lassoData.res.perfEnc],size(chanceData.perfEncAll,1),1) > chanceData.perfEncAll);
    auc_pEnc               = mean(repmat([lassoData.res.AUC_enc],size(chanceData.aucEncAll,1),1) > chanceData.aucEncAll);
    aucBest_pEnc           = mean(repmat([lassoData.resBest.AUC_enc],size(chanceData.aucBestEnc_all,1),1) > chanceData.aucBestEnc_all);
    perfEnc_p_all(s,:)     = perf_pEnc;
    aucEnc_p_all(s,:)      = auc_pEnc;
    aucBestEnc_p_all(s,:)  = aucBest_pEnc;
    perfEnc_all(s,:)       = [lassoData.res.perfEnc];
    aucEnc_all(s,:)        = [lassoData.res.AUC_enc];
    aucBestEnc_all(s,:)    = [lassoData.resBest.AUC_enc];
    lassoData.perfEnc      = [lassoData.res.perfEnc];
    lassoData.AUC_enc      = [lassoData.res.AUC_enc];
            
    % recalled vs not recalled vector
    rec       = vertcat(lassoData.Y);
    for t = 1:length(lassoData.res)
        
        
        %% YC1 direct time bin mapping
        classProb = lassoData.res(t).yPred;
        [~,ind] = sort(classProb);
        recSort = rec(ind);
        
        % now bin the sorted recall vector
        start = 1:(length(recSort)/3):length(recSort);
        stop = [start(2:end)-1 length(recSort)];        
        bins = NaN(1,length(stop));
        for r = 1:length(stop)
            bins(r) = mean(recSort(start(r):stop(r)));
        end
        quarts_all(t,:,s) = bins;
        
        %% YC to encoding mapping
        classProb = lassoData.res(t).yPredEnc;                
        [~,ind] = sort(classProb);
        recSort = rec(ind);
        
        % now bin the sorted recall vector
        bins = NaN(1,length(stop));
        for r = 1:length(stop)
            bins(r) = mean(recSort(start(r):stop(r)));
        end
        quartsEnc_all(t,:,s) = bins;
        
        %% YC1 best time bin mapping
        classProb = lassoData.resBest(t).yPred;
        [~,ind] = sort(classProb);
        recSort = rec(ind);
        
        % now bin the sorted recall vector
        bins = NaN(1,length(stop));
        for r = 1:length(stop)
            bins(r) = mean(recSort(start(r):stop(r)));
        end
        quartsBest_all(t,:,s) = bins;    
        
        %% YC1 best time bin mapping to encoding. I realize each timebin is the same...
        classProb = lassoData.resBest(t).yPredEnc;
        [~,ind] = sort(classProb);
        recSort = rec(ind);
        
        % now bin the sorted recall vector
        bins = NaN(1,length(stop));
        for r = 1:length(stop)
            bins(r) = mean(recSort(start(r):stop(r)));
        end
        quartsBestEnc_all(t,:,s) = bins;            
        
    end
    
    %% FIGURE 1a and b - classifier accuracy and AUC over time
    
    % create x labels based on time bins
    xBins    = round((params.timeBins) / 10) / 100;
    xBinsStr = {};
    for x = 1:size(xBins,1)
        xBinsStr{x} = [num2str(xBins(x,1)), '-', num2str(xBins(x,2))];
    end
    
    ylabels   = {'Classifier Accuracy (%)','Classifier AUC','Classifier AUC'};
    ps        = {[1-perf_p;1-perf_pEnc],[1-auc_p;1-auc_pEnc],[1-aucBest_p;1-aucBest_pEnc]};    
    fields    = {{'perf','perfEnc'},{'AUC','AUC_enc'},{'AUC_best','AUC_bestEnc'}};
    axPos     = [.15 .55 .8 .4;.15 .15 .8 .4];
    
    for i = 1:3
        fname = fullfile(figDir,[subj '_' fields{i}{1} '.eps']);
        figs_subj.(fields{i}{1}) = fname;
        if exist(fname,'file') && ~overwrite
            continue
        end
        figure(1)
        clf
        
        for j = 1:2
            figure(1)
            % first plot all the points as black
            axes('position',axPos(j,:));
            plot(lassoData.(fields{i}{j}),'k.','markersize',30,'linewidth',2.5)
            hold on
            
            % plot any significant time points as red. Bigger red dot indicates
            % the timepoint survived bonferroni correction. Small red dot is p
            % < .05
            p = ps{i}(j,:);
            thresh = .05/length(lassoData.(fields{i}{j}));
            h=plot(find(p < .05),lassoData.(fields{i}{j})(p < .05),'r.','markersize',30);
            set(h,'color',[140 15 15]/255)
            plot(find(p < thresh),lassoData.(fields{i}{j})(p < thresh),'r.','markersize',55)
            set(h,'color',[200 100 100]/255)
            
            % set axis and labels
            set(gca,'xtick',1:length(lassoData.(fields{i}{j})));
            set(gca,'xlim',[0 length(lassoData.(fields{i}{j}))+1]);
            set(gca,'xticklabel',xBinsStr);
            set(gca,'ylim',[0 1]);
            set(gca,'ytick',.2:.2:.8);
            if i == 1
                set(gca,'yticklabel',20:20:80)
            end
            
            if j == 2
                xlabel('Time (s)','fontsize',16)
                h=text(-1.3,1,ylabels{i},'fontsize',16);
                len = get(h,'extent');
                len = len(2);
                pos = get(h,'position');
                set(h,'position',[pos(1) pos(2)-len*(.5) 0])
                set(h,'rotation',90)
            end
            grid on
            set(gca,'gridlinestyle',':');
            set(gca,'fontsize',16)
            
            if j == 1
                fillX = [.5 .5 length(lassoData.(fields{i}{j}))+.05 length(lassoData.(fields{i}{j}))+.05];
                h=fill(fillX,[.91 .98 .98 .91],'w');
                set(h,'linestyle','none');
                
                % label each time point with the behavior
                labels = params.timeBinLabels;
                if isempty(labels);labels=repmat({''},1,size(params.timeBins,1));end
                for t = 1:length(lassoData.(fields{i}{j}))
                    start = t - .25;
                    stop  = t + .25;
                    plot([start stop],[.9 .9],'-k','linewidth',2)
                    h=text(t,.95,labels{t},'fontsize',12);
                    len = get(h,'extent');
                    len = len(3);
                    pos = get(h,'position');
                    set(h,'position',[pos(1)-len/2+.03 pos(2) 0]);
                end
            end
            
            % plot dashed line at 50%
            xlim = get(gca,'xlim');
            plot(xlim,[.5 .5],'--k','linewidth',1.5)
            box on
            
        end
        print('-depsc2','-tiff','-loose',fname);
        
        
        %% FIGURE 2 a and b - quartile plot for most signficant time period for perf and AUC
        fname = fullfile(figDir,[subj '_' fields{i}{1} '_quart.eps']);
        figs_subj.(['quarts',fields{i}{1}]) = fname;
        figure(2)
        clf        
        
        for j = 1:2
            figure(2)
            % compute quartiles
            [minP,~] = min(ps{i}(j,:));
            allMinPs = find(ps{i}(j,:)==minP);            
            [~,ind] = max(lassoData.(fields{i}{j})(allMinPs));
            bestTime = allMinPs(ind);
            
            %             best_time(s) = bestTime;
            for t = 1:length(lassoData.res)                
                if t == bestTime
                    quarts = quarts_all(t,:,s);
                end
            end
            
            if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
                axes('position',[axPos(j,1:3) .3]);
                hold on
                % plot quartiles based on change from mean
                bar(((quarts -(mean(rec)))/mean(rec))*100,'w','linewidth',2);
                %                 ylabel('Recall Change (%)','fontsize',16)
                set(gca,'fontsize',16)
                set(gca,'ylim',[-50 50])
                set(gca,'xlim',[0 5])
                set(gca,'xtick',1:4)
                grid on
                set(gca,'gridlinestyle',':');
                hold on
                if j == 2
                    xlabel('Quartile of Classifier Estimate','fontsize',16)
                    h=text(-.65,50,'Recall Change (%)','fontsize',16);
                    len = get(h,'extent');
                    len = len(2);
                    pos = get(h,'position');
                    set(h,'position',[pos(1) pos(2)-len/.65 0]);
                    set(h,'rotation',90)
                else
                    set(gca,'xticklabel','')
                end
                h=title([labels{bestTime} ' Period'],'fontsize',16);
                set(h,'fontweight','normal');                
            end
        end
        print('-depsc2','-tiff','-loose',fname);
    end
   
    figs = [figs;figs_subj];
    
end


% also make group plots/report
fprintf('Creating group plots.\n');
figs_group = [];
figs_group.quarts   = {};
figs_group.quarts_enc = {};
figs_group.quartsBest = {};
figs_group.quartsBestEnc = {};
figs_group.auc_hist = {};
figs_group.aucEnc_hist = {};
figs_group.aucBestEnc_hist = {};
figs_group.aucBest_hist = {};
figs_group.quartsBest = {};
figs_group.quartsBestEnc = {};

% compute average quartile measure for each time bin
meanRec_subj    = repmat(nanmean(quarts_all,2),[1,3,1]);
nSubj           = sum(~isnan(quarts_all),3);
figs_group.N    = nSubj(:,1);
quarts_err      = nanstd((quarts_all - meanRec_subj)./meanRec_subj,[],3)./sqrt(nSubj-1);
quarts_group    = nanmean((quarts_all - meanRec_subj)./meanRec_subj,3);
quartsEnc_err   = nanstd((quartsEnc_all - meanRec_subj)./meanRec_subj,[],3)./sqrt(nSubj-1);
quartsEnc_group = nanmean((quartsEnc_all - meanRec_subj)./meanRec_subj,3);

quartsBest_all_err   = nanstd((quartsBest_all - meanRec_subj)./meanRec_subj,[],3)./sqrt(nSubj-1);
quartsBest_all_group = nanmean((quartsBest_all - meanRec_subj)./meanRec_subj,3);
quartsBestEnc_all_err   = nanstd((quartsBestEnc_all - meanRec_subj)./meanRec_subj,[],3)./sqrt(nSubj-1);
quartsBestEnc_all_group = nanmean((quartsBestEnc_all - meanRec_subj)./meanRec_subj,3);


% labels for plotting
labels = params.timeBinLabels;
if isempty(labels);labels=repmat({''},1,size(quarts_group,1));end

% plot each time bin separately
for t = 1:size(quarts_group,1)
    
    %% QUARTILE PLOT - YC1 direct time bin mapping
    fname = fullfile(figDir,['group_quart_' labels{t} '.eps']);
    figs_group.quarts{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(2)
        clf
        bar(quarts_group(t,:)*100,'w','linewidth',2);
        hold on
        errorbar(1:3,quarts_group(t,:)*100,quarts_err(t,:)*196,'k','linewidth',2,'linestyle','none')
        
        xlabel('Quartile of Classifier Estimate','fontsize',20)
        ylabel('Recall Change (%)','fontsize',20)
        set(gca,'fontsize',20)
        set(gca,'ylim',[-100 100])
        set(gca,'xlim',[0 5])
        set(gca,'ylim',[-25 25])
        grid on
        set(gca,'gridlinestyle',':');
        hold on
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');
        print('-depsc2','-tiff','-loose',fname);
    end
    
    %% QUARTILE PLOT - YC1 time bin to encoding mapping
    fname = fullfile(figDir,['group_quart_enc_' labels{t} '.eps']);
    figs_group.quarts_enc{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(2)
        clf
        bar(quartsEnc_group(t,:)*100,'w','linewidth',2);
        hold on
        errorbar(1:3,quartsEnc_group(t,:)*100,quartsEnc_err(t,:)*196,'k','linewidth',2,'linestyle','none')
        
        xlabel('Quartile of Classifier Estimate','fontsize',20)
        ylabel('Recall Change (%)','fontsize',20)
        set(gca,'fontsize',20)
        set(gca,'ylim',[-100 100])
        set(gca,'xlim',[0 5])
        set(gca,'ylim',[-25 25])
        grid on
        set(gca,'gridlinestyle',':');
        hold on
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');
        print('-depsc2','-tiff','-loose',fname);
    end    
    
   %% QUARTILE PLOT - YC1 best to each mapping
    fname = fullfile(figDir,['group_quartBest_' labels{t} '.eps']);
    figs_group.quartsBest{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(2)
        clf
        bar(quartsBest_all_group(t,:)*100,'w','linewidth',2);
        hold on
        errorbar(1:3,quartsBest_all_group(t,:)*100,quartsBest_all_err(t,:)*196,'k','linewidth',2,'linestyle','none')
        
        xlabel('Quartile of Classifier Estimate','fontsize',20)
        ylabel('Recall Change (%)','fontsize',20)
        set(gca,'fontsize',20)
        set(gca,'ylim',[-100 100])
        set(gca,'xlim',[0 5])
        set(gca,'ylim',[-25 25])
        grid on
        set(gca,'gridlinestyle',':');
        hold on
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');
        print('-depsc2','-tiff','-loose',fname);
    end    
    
      %% QUARTILE PLOT - YC1 best to Enc mapping. 
    fname = fullfile(figDir,['group_quartBestEnc_' labels{t} '.eps']);
    figs_group.quartsBestEnc{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(2)
        clf
        bar(quartsBestEnc_all_group(t,:)*100,'w','linewidth',2);
        hold on
        errorbar(1:3,quartsBestEnc_all_group(t,:)*100,quartsBestEnc_all_err(t,:)*196,'k','linewidth',2,'linestyle','none')
        
        xlabel('Quartile of Classifier Estimate','fontsize',20)
        ylabel('Recall Change (%)','fontsize',20)
        set(gca,'fontsize',20)
        set(gca,'ylim',[-100 100])
        set(gca,'xlim',[0 5])
        set(gca,'ylim',[-25 25])
        grid on
        set(gca,'gridlinestyle',':');
        hold on
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');
        print('-depsc2','-tiff','-loose',fname);
    end     
    
    
    %% AUC HISTOGRAM PLOT DIRECT TIME BIN MAPPING
    fname = fullfile(figDir,['auc_hist_' labels{t} '.eps']);
    figs_group.auc_hist{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(3)
        clf
        auc = auc_all(:,t);
        sig  = auc_p_all(:,t) > .95;
        n1   = histc(auc(sig),0.025:.05:.975);
        % wtf
        if isrow(n1);n1=n1';end        
        n2   = histc(auc(~sig),0.025:.05:.975);
        h    = bar([.05:.05:1]*100,[n1 n2],1,'stacked','linewidth',2);
        xlabel('Classifier AUC','Fontsize',20);
        set(gca,'xlim',[.2 .8]*100);
        set(gca,'xlim',[0 100]);
        set(gca,'xtick',0:25:100)
        set(gca,'ylim',[0 15]);
        ylabel('Subject Count','Fontsize',20)
        set(h(2),'FaceColor','w');
        set(h(1),'FaceColor',[226 55 67]/255);
        grid on
        set(gca,'fontsize',20)
        set(gca,'gridlinestyle',':');
        box on
        hold on
        plot([50 50],[0 15],'--k','linewidth',2)
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');     
        print('-depsc2','-tiff','-loose',fname);              
    end     
    
    %% AUC HISTOGRAM PLOT TIME BIN TO 0-5 ENCODING
    fname = fullfile(figDir,['auc_histEnc_' labels{t} '.eps']);
    figs_group.aucEnc_hist{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(3)
        clf
        auc = aucEnc_all(:,t);
        sig  = aucEnc_p_all(:,t) > .95;
        n1   = histc(auc(sig),0.025:.05:.975);
        % wtf
        if isrow(n1);n1=n1';end
        n2   = histc(auc(~sig),0.025:.05:.975);
        h    = bar([.05:.05:1]*100,[n1 n2],1,'stacked','linewidth',2);
        xlabel('Classifier AUC','Fontsize',20);
        set(gca,'xlim',[.2 .8]*100);
        set(gca,'xlim',[0 100]);
        set(gca,'xtick',0:25:100)
        set(gca,'ylim',[0 15]);
        ylabel('Subject Count','Fontsize',20)
        set(h(2),'FaceColor','w');
        set(h(1),'FaceColor',[226 55 67]/255);
        grid on
        set(gca,'fontsize',20)
        set(gca,'gridlinestyle',':');
        box on
        hold on
        plot([50 50],[0 15],'--k','linewidth',2)
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');     
        print('-depsc2','-tiff','-loose',fname);        
    end        
    
    %% AUC HISTOGRAM PLOT BEST TO EACH
    fname = fullfile(figDir,['aucBest_hist_' labels{t} '.eps']);
    figs_group.aucBest_hist{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(3)
        clf
        auc = aucBest_all(:,t);
        sig  = aucBest_p_all(:,t) > .95;
        n1   = histc(auc(sig),0.025:.05:.975);
        % wtf
        if isrow(n1);n1=n1';end        
        n2   = histc(auc(~sig),0.025:.05:.975);
        h    = bar([.05:.05:1]*100,[n1 n2],1,'stacked','linewidth',2);
        xlabel('Classifier AUC','Fontsize',20);
        set(gca,'xlim',[.2 .8]*100);
        set(gca,'xlim',[0 100]);
        set(gca,'xtick',0:25:100)
        set(gca,'ylim',[0 15]);
        ylabel('Subject Count','Fontsize',20)
        set(h(2),'FaceColor','w');
        set(h(1),'FaceColor',[226 55 67]/255);
        grid on
        set(gca,'fontsize',20)
        set(gca,'gridlinestyle',':');
        box on
        hold on
        plot([50 50],[0 15],'--k','linewidth',2)
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');     
        print('-depsc2','-tiff','-loose',fname);              
    end        
    
       %% AUC HISTOGRAM PLOT BEST TO ENC
    fname = fullfile(figDir,['aucBestEnc_hist_' labels{t} '.eps']);
    figs_group.aucBestEnc_hist{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(3)
        clf
        auc = aucBestEnc_all(:,t);
        sig  = aucBestEnc_p_all(:,t) > .95;
        n1   = histc(auc(sig),0.025:.05:.975);
        % wtf
        if isrow(n1);n1=n1';end        
        n2   = histc(auc(~sig),0.025:.05:.975);
        h    = bar([.05:.05:1]*100,[n1 n2],1,'stacked','linewidth',2);
        xlabel('Classifier AUC','Fontsize',20);
        set(gca,'xlim',[.2 .8]*100);
        set(gca,'xlim',[0 100]);
        set(gca,'xtick',0:25:100)
        set(gca,'ylim',[0 15]);
        ylabel('Subject Count','Fontsize',20)
        set(h(2),'FaceColor','w');
        set(h(1),'FaceColor',[226 55 67]/255);
        grid on
        set(gca,'fontsize',20)
        set(gca,'gridlinestyle',':');
        box on
        hold on
        plot([50 50],[0 15],'--k','linewidth',2)
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');     
        print('-depsc2','-tiff','-loose',fname);              
    end    
    
end

% AUC OVER DIRECT TIME BIN MAPPING
fname = fullfile(figDir,'auc_time_direct.eps');
figs_group.auc_time_direct = fname;
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(3)
    clf
    [h,p] = ttest(auc_all,.5);
    sigCorr = p*size(auc_all,2) < .05;
    sig     = p <.05 & ~sigCorr;
    h=bar(find(~sig & ~sigCorr),nanmean(auc_all(:,~sig & ~sigCorr)),'w','linewidth',2);
    set(h,'facecolor',[.5 .5 .5])
    hold on
    if any(sig)
        h=bar(find(sig),nanmean(auc_all(:,sig)),'w','linewidth',2);
        set(gca,'ylim',[.4 .6])
        set(h,'facecolor',[200 100 100]/255)
    end
    if any(sigCorr)
        h=bar(find(sigCorr),nanmean(auc_all(:,sigCorr)),'w','linewidth',2);
        set(gca,'ylim',[.4 .6])
        set(h,'facecolor',[140 15 15]/255)
    end
    err = nanstd(auc_all)./sqrt(sum(~isnan(auc_all))-1);
    errorbar(1:size(auc_all,2),nanmean(auc_all),err*1.96,'k','linewidth',2,'linestyle','none')
    plot([0 size(auc_all,2)+1],[.5 .5],'--k','linewidth',2)
    grid on
    set(gca,'gridlinestyle',':');
    set(gca,'xlim',[0 size(auc_all,2)+1]);
    set(gca,'xtick',1:size(auc_all,2));
    set(gca,'xticklabel',xBinsStr);
    set(gca,'ylim',[.4 .6])
    set(gca,'ytick',.4:.05:.6)
    ylabel('Classifier AUC','Fontsize',16);
    xlabel('Timebin','Fontsize',16);
    set(gca,'fontsize',16)
    print('-depsc2','-tiff','-loose',fname);   
end

% AUC OVER TIME BIN TO ALL ENC MAPPING
fname = fullfile(figDir,'auc_time_indirect.eps');
figs_group.auc_time_indirect = fname;
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(4)
    clf
    [h,p] = ttest(aucEnc_all,.5);
    sigCorr = p*size(aucEnc_all,2) < .05;
    sig     = p <.05 & ~sigCorr;
    h=bar(find(~sig & ~sigCorr),nanmean(aucEnc_all(:,~sig & ~sigCorr)),'w','linewidth',2);
    set(h,'facecolor',[.5 .5 .5])
    hold on
    if any(sig)
        h=bar(find(sig),nanmean(aucEnc_all(:,sig)),'w','linewidth',2);
        set(gca,'ylim',[.4 .6])
        set(h,'facecolor',[200 100 100]/255)
    end
    if any(sigCorr)
        h=bar(find(sigCorr),nanmean(aucEnc_all(:,sigCorr)),'w','linewidth',2);
        set(gca,'ylim',[.4 .6])
        set(h,'facecolor',[140 15 15]/255)
    end
    err = nanstd(aucEnc_all)./sqrt(sum(~isnan(aucEnc_all))-1);
    errorbar(1:size(aucEnc_all,2),nanmean(aucEnc_all),err*1.96,'k','linewidth',2,'linestyle','none')
    plot([0 size(aucEnc_all,2)+1],[.5 .5],'--k','linewidth',2)
    grid on
    set(gca,'gridlinestyle',':');
    set(gca,'xlim',[0 size(aucEnc_all,2)+1]);
    set(gca,'xtick',1:size(aucEnc_all,2));
    set(gca,'xticklabel',xBinsStr);
    set(gca,'ylim',[.4 .6])
    set(gca,'ytick',.4:.05:.6)
    ylabel('Classifier AUC','Fontsize',16);
    xlabel('Timebin','Fontsize',16);
    set(gca,'fontsize',16)
    print('-depsc2','-tiff','-loose',fname);   
end

% AUC BEST TO EACH
fname = fullfile(figDir,'aucBest_time_direct.eps');
figs_group.aucBest_time_direct = fname;
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(3)
    clf
    [h,p] = ttest(aucBest_all,.5);
    sigCorr = p*size(aucBest_all,2) < .05;
    sig     = p <.05 & ~sigCorr;
    h=bar(find(~sig & ~sigCorr),nanmean(aucBest_all(:,~sig & ~sigCorr)),'w','linewidth',2);
    set(h,'facecolor',[.5 .5 .5])
    hold on
    if any(sig)
        h=bar(find(sig),nanmean(aucBest_all(:,sig)),'w','linewidth',2);
        set(gca,'ylim',[.4 .6])
        set(h,'facecolor',[200 100 100]/255)
    end
    if any(sigCorr)
        h=bar(find(sigCorr),nanmean(aucBest_all(:,sigCorr)),'w','linewidth',2);
        set(gca,'ylim',[.4 .6])
        set(h,'facecolor',[140 15 15]/255)
    end
    err = nanstd(aucBest_all)./sqrt(sum(~isnan(aucBest_all))-1);
    errorbar(1:size(aucBest_all,2),nanmean(aucBest_all),err*1.96,'k','linewidth',2,'linestyle','none')
    plot([0 size(aucBest_all,2)+1],[.5 .5],'--k','linewidth',2)
    grid on
    set(gca,'gridlinestyle',':');
    set(gca,'xlim',[0 size(aucBest_all,2)+1]);
    set(gca,'xtick',1:size(aucBest_all,2));
    set(gca,'xticklabel',xBinsStr);
    set(gca,'ylim',[.4 .6])
    set(gca,'ytick',.4:.05:.6)
    ylabel('Classifier AUC','Fontsize',16);
    xlabel('Timebin','Fontsize',16);
    set(gca,'fontsize',16)
    print('-depsc2','-tiff','-loose',fname);   
end

% AUC BEST TO ENC
fname = fullfile(figDir,'aucBestEnc_time_direct.eps');
figs_group.aucBestEnc_time_direct = fname;
if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
    figure(3)
    clf
    [h,p] = ttest(aucBestEnc_all,.5);
    sigCorr = p*size(aucBestEnc_all,2) < .05;
    sig     = p <.05 & ~sigCorr;
    h=bar(find(~sig & ~sigCorr),nanmean(aucBestEnc_all(:,~sig & ~sigCorr)),'w','linewidth',2);
    set(h,'facecolor',[.5 .5 .5])
    hold on
    if any(sig)
        h=bar(find(sig),nanmean(aucBestEnc_all(:,sig)),'w','linewidth',2);
        set(gca,'ylim',[.4 .6])
        set(h,'facecolor',[200 100 100]/255)
    end
    if any(sigCorr)
        h=bar(find(sigCorr),nanmean(aucBestEnc_all(:,sigCorr)),'w','linewidth',2);
        set(gca,'ylim',[.4 .6])
        set(h,'facecolor',[140 15 15]/255)
    end
    err = nanstd(aucBestEnc_all)./sqrt(sum(~isnan(aucBestEnc_all))-1);
    errorbar(1:size(aucBestEnc_all,2),nanmean(aucBestEnc_all),err*1.96,'k','linewidth',2,'linestyle','none')
    plot([0 size(aucBestEnc_all,2)+1],[.5 .5],'--k','linewidth',2)
    grid on
    set(gca,'gridlinestyle',':');
    set(gca,'xlim',[0 size(aucBestEnc_all,2)+1]);
    set(gca,'xtick',1:size(aucBestEnc_all,2));
    set(gca,'xticklabel',xBinsStr);
    set(gca,'ylim',[.4 .6])
    set(gca,'ytick',.4:.05:.6)
    ylabel('Classifier AUC','Fontsize',16);
    xlabel('Timebin','Fontsize',16);
    set(gca,'fontsize',16)
    print('-depsc2','-tiff','-loose',fname);   
end

keyboard
good = ~cellfun('isempty',{figs.subj});
figs = figs(good);
texName = 'YC2_lassoChance_report.tex';
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
texName = 'YC2_group_lassoChance_report.tex';
write_texfile_group(saveDir,texName,figs_group)

curr_dir = pwd;
cd(saveDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(saveDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);




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
    fprintf(fid,'\\includegraphics[width=0.45\\textwidth]{%s}\n',figs(s).perf);
    fprintf(fid,'\\includegraphics[width=0.45\\textwidth]{%s}\n',figs(s).AUC);
%     fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).quartsperf);
%     fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).quartsAUC);
    fprintf(fid,'\\caption{%s: region: %s, %d electrodes}\n\n',strrep(figs(s).subj,'_',' '),figs(s).region,figs(s).nElecs);
    fprintf(fid,'\\end{figure}\n\n\n');
    if mod(s,3) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
end 

fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);
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
fprintf(fid,'\\rhead{YC1 Group Report Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.quarts{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject average quartile by time bin. YC1 model from each time bin is applied to the \\textbf{same} time bin in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');





fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.quarts_enc{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject average quartile by time bin. YC1 model from each time bin is applied to the \\textbf{0-5 encoding bin} in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');


fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.quartsBest{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject average quartile by time bin. YC1 model from BEST time bin is applied to the \\textbf{same} time bin in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');

fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.quartsBestEnc{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject average quartile by time bin. YC1 model from BEST time bin is applied to the \\textbf{ENC} time bin in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');


fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.auc_hist{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject AUC histogram by time bin. YC1 model from each time bin is applied to the \\textbf{same} time bin in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');


fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.aucEnc_hist{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject AUC histogram by time bin. YC1 model from each time bin is applied to the \\textbf{0-5 encoding bin} in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');



fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.aucBest_hist{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject AUC histogram by time bin. YC1 model from BEST time bin is applied to the \\textbf{same} time bin in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');


fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.aucBestEnc_hist{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject AUC histogram by time bin. YC1 model from BEST time bin is applied to the \\textbf{0-5 encoding bin} in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');



fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');

fprintf(fid,'\\includegraphics[width=0.45\\textwidth]{%s}\n',figs.auc_time_direct);
fprintf(fid,'\\includegraphics[width=0.45\\textwidth]{%s}\n',figs.auc_time_indirect);
fprintf(fid,'\\caption{%d Subjects: Average AUC over time. Dark red: significant after correcting for number of time bins. Light red: $p<.05$. \\textbf{Left:} YC1 model for each time bin applied to the same YC2 timebin. \\textbf{Right:} YC1 model for each time bin applied to average 0-5 second bin.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');

fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');

fprintf(fid,'\\includegraphics[width=0.45\\textwidth]{%s}\n',figs.aucBest_time_direct);
fprintf(fid,'\\includegraphics[width=0.45\\textwidth]{%s}\n',figs.aucBestEnc_time_direct);
fprintf(fid,'\\caption{%d Subjects: Average AUC over time. Dark red: significant after correcting for number of time bins. Light red: $p<.05$. \\textbf{Left:} YC1 model for each time bin applied to the same YC2 timebin. \\textbf{Right:} YC1 model for each time bin applied to average 0-5 second bin.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');

fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);






