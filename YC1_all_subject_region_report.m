function YC1_all_subject_region_report()


% analysis settings
% -----------------

% contrast to do
anas = {'correct_incorrect_resids'};
ana_names = {'Correct - Incorrect'};

% limit to just correct trials?
only_correct = false;

% save directory
saveDirs = cell(1,length(anas));
prefix = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TOM: Path is set here %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for anNum = 1:length(anas)
    saveDirs{anNum} = fullfile(prefix,'/data10/scratch/jfm2/YC1/group',anas{anNum});
    if ~exist(saveDirs{anNum},'dir')
        error('Invalid analysis: %s',anas{anNum})
    elseif ~exist(fullfile(saveDirs{anNum},'figs'),'dir')
        mkdir(fullfile(saveDirs{anNum},'figs'));
    end
end

% get list of YC subjects
subjs = get_subs('RAM_YC1');

% % subjs = subjs(~strcmp(subjs,'TJ078'));
process_subj(subjs,anas,saveDirs,ana_names);


function process_subj(subj,anas,saveDirs,ana_names)

% load task config to get the frequencies
config = RAM_config('RAM_YC1');

% Pick an roi 'hipp','ec','mtl','frontal','parietal','temporal', 'occipital','limbic'
rois = {'hipp'};
figs = struct([]);
if isempty(strfind(version(),'2014'));
    fixFiguresFor2014 = false;
    lineWidth = 2;
else
    fixFiguresFor2014 = true;
    lineWidth = 4;
end

for anNum = 1:length(anas)
    
    saveDir = saveDirs{anNum};
    ana = anas{anNum};
    figs(anNum).ana_name = ana_names{anNum};
    [~,ana_path] = fileparts(saveDir);
    figs(anNum).ana_path = ana_path;
    
    % will hold combined results
    sig_countsLTA = NaN(length(rois),2);
    sig_countsHTA = NaN(length(rois),2);
    tstat_LTA = NaN(length(rois),2);
    tstat_HTA = NaN(length(rois),2);
    corr_LTA = NaN(length(rois),2);
    corr_HTA = NaN(length(rois),2);
    corr_G = NaN(length(rois),2);
    corr_HFA = NaN(length(rois),2);
    sig_countsG = NaN(length(rois),2);
    sig_countsHFA = NaN(length(rois),2);
    tstat_G = NaN(length(rois),2);
    tstat_HFA = NaN(length(rois),2);
    num_elecs = NaN(length(rois),1);
    labels = {};
    link_text = {};
    
    if ~fixFiguresFor2014
        figure('units','normalized','paperpositionmode','auto','position',[0.7504    0.6850    0.2188    0.2662])
    else
        figure('units','normalized','paperpositionmode','auto','outerposition',[0.7504    0.6850    0.2188    0.2])
    end
    for r = 1:length(rois)
        
        disp(rois{r})
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% THIS WHOLE MESS LOADS ALL THE DATA CREATED BY        %%%
        %%% YC1_SUBJECT_SUMMARY AND COMBINES INTO ONE STRUCTURE  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        all_subjs = [];
        for s = 1:length(subj)
            if exist(fullfile(saveDir,[subj{s},'_',rois{r},'_pow.mat']),'file')
                regionData = load(fullfile(saveDir,[subj{s},'_',rois{r},'_pow.mat']));
                regionData.er = repmat(regionData.er,[size(regionData.powLTA,1),1]);
                
                % this is hack for now because the number of trials is
                % not consistent for across subjects. assuming no
                % patients has more than 1000 trals
                powLTA = NaN(size(regionData.powLTA,1),1000);
                powLTA(:,1:size(regionData.powLTA,2)) = regionData.powLTA;
                regionData.powLTA = powLTA;
                
                powHTA = NaN(size(regionData.powHTA,1),1000);
                powHTA(:,1:size(regionData.powHTA,2)) = regionData.powHTA;
                regionData.powHTA = powHTA;
                
                powHFA = NaN(size(regionData.powHTA,1),1000);
                powHFA(:,1:size(regionData.powHTA,2)) = regionData.powHTA;
                regionData.powHFA = powHFA;
                
                powG = NaN(size(regionData.powHTA,1),1000);
                powG(:,1:size(regionData.powHTA,2)) = regionData.powHTA;
                regionData.powG = powG;
                
                er = NaN(size(regionData.er,1),1000);
                er(:,1:size(regionData.er,2)) = regionData.er;
                regionData.er = er;
                
                if isempty(all_subjs)
                    fnames = {'powCond1ByElec','powCond2ByElec',...
                        'statsG','statsHFA','statsLTA','statsHTA',...
                        'rLTA','pLTA','rHTA','pHTA','powLTA',...
                        'powHTA','er','rG','pG','rHFA','pHFA','powHFA'};
                    all_subjs = cell2struct(cell(1,length(fnames)),{fnames{:}},2);
                end
                
                for f = 1:length(fnames)
                    if iscolumn(regionData.(fnames{f})) && length(regionData.statsLTA) > 1
                        if ~ismember(fnames{f},{'powHTA','powLTA','er'})
                            regionData.(fnames{f}) = regionData.(fnames{f})';
                        end
                    end
                    if ismember(fnames{f},{'powHTA','powLTA','er','powG','powHFA'})
                        all_subjs.(fnames{f}) = cat(1,all_subjs.(fnames{f}),regionData.(fnames{f}));
                    else
                        all_subjs.(fnames{f}) = cat(2,all_subjs.(fnames{f}),regionData.(fnames{f}));
                    end
                end
            end
            
        end
        regionData = all_subjs;
        
        
        num_elecs(r) = length(regionData.statsLTA);
        labels{r} = sprintf('%s (%d)',rois{r},num_elecs(r));
        link_text{end+1} = rois{r};
        
        % low theta data
        ltaSig = [regionData.statsLTA.p] < .05;
        ltaCountPos = sum([regionData.statsLTA(ltaSig).tstat] > 0);
        sig_countsLTA(r,1) = ltaCountPos/num_elecs(r);
        sig_countsLTA(r,2) = (sum(ltaSig) - ltaCountPos)/num_elecs(r);
        tstat_LTA(r,1) = nanmean([regionData.statsLTA.tstat]);
        tstat_LTA(r,2) = 1 * nanstd([regionData.statsLTA.tstat])/(sqrt(num_elecs(r))-1);
        corr_LTA(r,1) = nanmean([regionData.rLTA]);
        corr_LTA(r,2) = 1 * nanstd([regionData.rLTA])/(sqrt(num_elecs(r))-1);
        
        % high theta data
        htaSig = [regionData.statsHTA.p] < .05;
        htaCountPos = sum([regionData.statsHTA(htaSig).tstat] > 0);
        sig_countsHTA(r,1) = htaCountPos/num_elecs(r);
        sig_countsHTA(r,2) = (sum(htaSig) - htaCountPos)/num_elecs(r);
        tstat_HTA(r,1) = nanmean([regionData.statsHTA.tstat]);
        tstat_HTA(r,2) = 1 * nanstd([regionData.statsHTA.tstat])/(sqrt(num_elecs(r))-1);
        corr_HTA(r,1) = nanmean([regionData.rHTA]);
        corr_HTA(r,2) = 1 * nanstd([regionData.rHTA])/(sqrt(num_elecs(r))-1);
        
        % gamma
        gSig = [regionData.statsG.p] < .05;
        gCountPos = sum([regionData.statsG(gSig).tstat] > 0);
        sig_countsG(r,1) = gCountPos/num_elecs(r);
        sig_countsG(r,2) = (sum(gSig) - gCountPos)/ ...
            num_elecs(r);
        tstat_G(r,1) = nanmean([regionData.statsG.tstat]);
        tstat_G(r,2) = 1 * nanstd([regionData.statsG.tstat])/(sqrt(num_elecs(r))-1);
        corr_G(r,1) = nanmean([regionData.rG]);
        corr_G(r,2) = 1 * nanstd([regionData.rG])/(sqrt(num_elecs(r))-1);
        
        % hfa
        hfaSig = [regionData.statsHFA.p] < .05;
        hfaCountPos = sum([regionData.statsHFA(hfaSig).tstat] > 0);
        sig_countsHFA(r,1) = hfaCountPos/num_elecs(r);
        sig_countsHFA(r,2) = (sum(hfaSig) - hfaCountPos)/ ...
            num_elecs(r);
        tstat_HFA(r,1) = nanmean([regionData.statsHFA.tstat]);
        tstat_HFA(r,2) = 1 * nanstd([regionData.statsHFA.tstat])/(sqrt(num_elecs(r))-1);
        corr_HFA(r,1) = nanmean([regionData.rHFA]);
        corr_HFA(r,2) = 1 * nanstd([regionData.rHFA])/(sqrt(num_elecs(r))-1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% also i'm plotting the difference between conditions at %%%
        %%%% every frequency                                        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [h,p,c,s] = ttest([regionData.powCond1ByElec - regionData.powCond2ByElec]');
        powerDiff = nanmean(regionData.powCond1ByElec - regionData.powCond2ByElec,2)';
        powerDiff_STD = nanstd(regionData.powCond1ByElec - regionData.powCond2ByElec,[],2)/sqrt(num_elecs(r)-1);
        
        plot([log10(config.distributedParams.freQ(1)) log10(config.distributedParams.freQ(end))],[0 0],'-k','linewidth',1)
        hold on
        plotErrorRegions(log10(config.distributedParams.freQ),powerDiff,powerDiff_STD,{[.5 .5 .5]});
        hold on
        
        % set pthresh for this plot
        pthresh = .05/length(config.distributedParams.freQ);
        pthresh = .01;
        if any(p < pthresh)
            colors = {[226,55,67]/255,[50,124,203]/255};
            sig = p < pthresh;
            if any(s.tstat > 0 & sig)
                plot(log10(config.distributedParams.freQ(s.tstat > 0 & sig)),.7,'.','markersize',6,'color',colors{1});
            end
            if any(s.tstat < 0 & sig)
                plot(log10(config.distributedParams.freQ(s.tstat < 0 & sig)),.7,'.','markersize',6,'color',colors{2});
            end
        end
        grid on
        
        set(gca,'xtick',log10(2.^(0:ceil(sqrt(max(config.distributedParams.freQ))))))
        set(gca,'xlim',[log10(config.distributedParams.freQ(1)) log10(config.distributedParams.freQ(end))])
        set(gca,'ylim',[-.5 .5])
        set(gca,'xticklabel','')
        title(labels{r})
        
        if fixFiguresFor2014
            set(gca,'titlefontweight','normal')
            set(gca,'GridLineStyle',':')
            set(gca,'GridColor','k')
            set(gca,'xticklabel',(2.^(0:ceil(sqrt(max(config.distributedParams.freQ))))))
            set(gca,'fontsize',7)
        end
        
        set(gca,'xticklabel',(2.^(0:ceil(sqrt(max(config.distributedParams.freQ))))))
        xlabel('Frequency (Hz)','fontsize',14)
        ylabel('Z(Power) Diff.','fontsize',14)
    end
    
    set(gcf,'PaperPositionMode','auto');
    figs(anNum).zpower = fullfile(saveDir,'figs',['all_subj_' 'zpow_' ana '.eps']);
    print(figs(anNum).zpower,'-depsc2');
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% SIG ELECTRODE COUNT BAR PLOT %%%
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %     if fixFiguresFor2014
    %         figure('units','normalized','paperpositionmode','auto','outerposition',[0.7504    0.6850    0.45    0.4],'position',[0.7504    0.6850    0.35    0.45])
    %     else
    %         figure(2)
    %     end
    %     clf
    %
    %     colors = {[226,55,67]/255,[50,124,203]/255};
    %     h = bar(1:3:22,sig_countsLTA*100,.25,'stacked','linewidth',lineWidth);
    %     set(h,{'FaceColor'},{colors{1};colors{2}});
    %     hold on
    %     h = bar(2:3:23,sig_countsHTA*100,.25,'stacked','linewidth',lineWidth);
    %     set(h,{'FaceColor'},{colors{1};colors{2}});
    %     set(gca,'xtick',1.5:3:23)
    %     set(gca,'xticklabel',labels)
    %     set(gca,'ylim',[0 110]);
    %     set(gca,'ytick',[0:20:100]);
    %
    %     ylabel('% Sig. Elecs.','fontsize',16)
    %     set(gca,'fontsize',16)
    %     grid on
    %     if fixFiguresFor2014
    %         set(gca,'XTickLabelRotation',-45)
    %         set(gca,'GridLineStyle',':')
    %         set(gca,'GridColor','k')
    %         ylabel('% Sig. Elecs.','fontsize',30)
    %         set(gca,'fontsize',30)
    %         set(gca,'clipping','off')
    %     else
    %         rotateXLabels(gca,45)
    %     end
    %
    %     set(gcf,'PaperPositionMode','auto');
    %     figs(anNum).sigCount_tTest = fullfile(saveDir,'figs',['all_subj_' 'tTest_pcount_' ana '.eps']);
    %     print(figs(anNum).sigCount_tTest,'-loose','-depsc2');
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% TSTAT BAR PLOT %%%
    %%%%%%%%%%%%%%%%%%%%%%
    if fixFiguresFor2014
        figure('units','normalized','paperpositionmode','auto','outerposition',[0.7504    0.6850    0.45    0.4],'position',[0.7504    0.6850    0.35    0.45])
    else
        figure(3)
    end
    clf
    hold on
    colors = {[226,55,67]/255,[50,124,203]/255};
    x1 = 1:5:(length(rois)*5);
    x2 = 2:5:(length(rois)*5+1);
    x3 = 3:5:(length(rois)*5+2);
    x4 = 4:5:(length(rois)*5+3);
    for i = 1:length(x1)
        c = colors{2};
        if tstat_LTA(i,1) > 0
            c = colors{1};
        end
        h=bar(x1(i),tstat_LTA(i,1),'linewidth',lineWidth);
        set(h,'FaceColor',c)
        errorbar(x1(i),tstat_LTA(i,1),tstat_LTA(i,2),'k','linewidth',lineWidth)
        
        c = colors{2};
        if tstat_HTA(i,1) > 0
            c = colors{1};
        end
        h=bar(x2(i),tstat_HTA(i,1),'linewidth',lineWidth);
        set(h,'FaceColor',c)
        errorbar(x2(i),tstat_HTA(i,1),tstat_HTA(i,2),'k','linewidth',lineWidth)
        
        c = colors{2};
        if tstat_G(i,1) > 0
            c = colors{1};
        end
        h=bar(x3(i),tstat_G(i,1),'linewidth',lineWidth);
        set(h,'FaceColor',c)
        errorbar(x3(i),tstat_G(i,1),tstat_G(i,2),'k','linewidth',lineWidth)
        
        c = colors{2};
        if tstat_HFA(i,1) > 0
            c = colors{1};
        end
        h=bar(x4(i),tstat_HFA(i,1),'linewidth',lineWidth);
        set(h,'FaceColor',c)
        errorbar(x4(i),tstat_HFA(i,1),tstat_HFA(i,2),'k','linewidth',lineWidth)
        
    end
    
    set(gca,'xtick',2.5:5:(length(rois)*5+2))
    set(gca,'xtick',1:4)
    %     set(gca,'xticklabel',labels)
    set(gca,'xticklabel',{'Low \theta','High \theta','\gamma','HFA'})
    
    ylabel('Mean t-stat','fontsize',16)
    set(gca,'fontsize',16)
    grid on
    box on
    if fixFiguresFor2014
        %         set(gca,'XTickLabelRotation',-45)
        set(gca,'GridLineStyle',':')
        set(gca,'GridColor','k')
        ylabel('Mean t-stat','fontsize',30)
        set(gca,'fontsize',30)
        set(gca,'clipping','off')
    else
        rotateXLabels(gca,45)
    end
    
    set(gcf,'PaperPositionMode','auto');
    figs(anNum).tTest_tStat = fullfile(saveDir,'figs',['all_subj_' 'tTest_tstat_' ana '.eps']);
    print(figs(anNum).tTest_tStat,'-loose','-depsc2');
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% CORR BAR PLOT %%%
    %%%%%%%%%%%%%%%%%%%%%%
    if fixFiguresFor2014
        figure('units','normalized','paperpositionmode','auto','outerposition',[0.7504    0.6850    0.45    0.4],'position',[0.7504    0.6850    0.35    0.45])
    else
        figure(4)
    end
    clf
    hold on
    colors = {[226,55,67]/255,[50,124,203]/255};
    colors = colors([2 1]);
    x1 = 1:5:(length(rois)*5);
    x2 = 2:5:(length(rois)*5+1);
    x3 = 3:5:(length(rois)*5+2);
    x4 = 4:5:(length(rois)*5+3);
    for i = 1:length(x1)
        c = colors{2};
        if corr_LTA(i,1) > 0
            c = colors{1};
        end
        h=bar(x1(i),corr_LTA(i,1),'linewidth',lineWidth);
        set(h,'FaceColor',c)
        errorbar(x1(i),corr_LTA(i,1),corr_LTA(i,2),'k','linewidth',lineWidth)
        
        c = colors{2};
        if corr_HTA(i,1) > 0
            c = colors{1};
        end
        h=bar(x2(i),corr_HTA(i,1),'linewidth',lineWidth);
        set(h,'FaceColor',c)
        errorbar(x2(i),corr_HTA(i,1),corr_HTA(i,2),'k','linewidth',lineWidth)
        
        c = colors{2};
        if corr_G(i,1) > 0
            c = colors{1};
        end
        h=bar(x3(i),corr_G(i,1),'linewidth',lineWidth);
        set(h,'FaceColor',c)
        errorbar(x3(i),corr_G(i,1),corr_G(i,2),'k','linewidth',lineWidth)
        
        c = colors{2};
        if corr_HFA(i,1) > 0
            c = colors{1};
        end
        h=bar(x4(i),corr_HFA(i,1),'linewidth',lineWidth);
        set(h,'FaceColor',c)
        errorbar(x4(i),corr_HFA(i,1),corr_HFA(i,2),'k','linewidth',lineWidth)
    end
    
    set(gca,'xtick',2.5:5:(length(rois)*5+2))
    set(gca,'xtick',1:4)
    %     set(gca,'xticklabel',labels)
    set(gca,'xticklabel',{'Low \theta','High \theta','\gamma','HFA'})
    
    ylabel('Mean Pearson Coef.','fontsize',16)
    set(gca,'fontsize',16)
    grid on
    box on
    if fixFiguresFor2014
        %         set(gca,'XTickLabelRotation',-45)
        set(gca,'GridLineStyle',':')
        set(gca,'GridColor','k')
        ylabel('Mean Pearson Coef.','fontsize',30)
        set(gca,'fontsize',30)
        set(gca,'clipping','off')
    else
        rotateXLabels(gca,45)
    end
    
    set(gcf,'PaperPositionMode','auto');
    figs(anNum).corrCoef = fullfile(saveDir,'figs',['all_subj_' 'corrCoef' ana '.eps']);
    print(figs(anNum).corrCoef,'-loose','-depsc2');
    keyboard
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SIG ELECTRODE CORR BAR PLOT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     if fixFiguresFor2014
    %         figure('units','normalized','paperpositionmode','auto','outerposition',[0.7504    0.6850   0.45    0.4],'position',[0.7504    0.6850    0.35    0.45])
    %     else
    %         figure(4)
    %     end
    %     clf
    %     colors = {[226,55,67]/255,[50,124,203]/255};
    %     h = bar(1:3:22,sig_corrCountsLTA*100,.25,'stacked','linewidth',lineWidth);
    %     set(h,{'FaceColor'},{colors{1};colors{2}});
    %     hold on
    %     h = bar(2:3:23,sig_corrCountsHTA*100,.25,'stacked','linewidth',lineWidth);
    %     set(h,{'FaceColor'},{colors{1};colors{2}});
    %     set(gca,'xtick',1.5:3:23)
    %     set(gca,'xticklabel',labels)
    %     set(gca,'ylim',[0 110]);
    %     set(gca,'ytick',[0:20:100]);
    %
    %     ylabel('% Sig. Elecs.','fontsize',16)
    %     set(gca,'fontsize',16)
    %     grid on
    %     if fixFiguresFor2014
    %         set(gca,'XTickLabelRotation',-45)
    %         set(gca,'GridLineStyle',':')
    %         set(gca,'GridColor','k')
    %         ylabel('% Sig. Elecs.','fontsize',30)
    %         set(gca,'fontsize',30)
    %         set(gca,'clipping','off')
    %     else
    %         rotateXLabels(gca,45)
    %     end
    %     set(gcf,'PaperPositionMode','auto');
    %     figs(anNum).sigCount_corr = fullfile(saveDir,'figs',['all_subj_' 'corr_' ana '.eps']);
    %     print(figs(anNum).sigCount_corr,'-loose','-depsc2');
    %
    %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%% PEAK FREQUENCY DIFF BAR PLOT %%%
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     if fixFiguresFor2014
    %         figure('units','normalized','paperpositionmode','auto','outerposition',[0.7504    0.6850   0.45    0.4],'position',[0.7504    0.6850    0.35    0.45])
    %     else
    %         figure(4)
    %     end
    %     clf
    %     bar(1:length(rois),freq_change,'facecolor','w','linewidth',lineWidth);
    %     hold on
    %     errorbar(1:length(rois),freq_change,freq_change_err,'k','linestyle','none','linewidth',lineWidth)
    %     % set(gca,'xtick',1:1:length(rois)+1)
    %     set(gca,'xticklabel',labels)
    %     % set(gca,'ylim',[0 110]);
    %     % set(gca,'ytick',[0:20:100]);
    %
    %     set(gca,'xlim',[.25 length(rois)+.75]);
    %     ylabel('Peak Freq. Diff. (Hz)','fontsize',16)
    %     grid on
    %     set(gca,'fontsize',16)
    %     if fixFiguresFor2014
    %         set(gca,'XTickLabelRotation',-45)
    %         set(gca,'GridLineStyle',':')
    %         set(gca,'GridColor','k')
    %         ylabel('Peak Freq. Diff. (Hz)','fontsize',30)
    %         set(gca,'fontsize',30)
    %         set(gca,'clipping','off')
    %     else
    %         rotateXLabels(gca,45)
    %     end
    %     set(gcf,'PaperPositionMode','auto');
    %     figs(anNum).peak_diff = fullfile(saveDir,'figs',['all_subj_' 'peak_diff_' ana '.eps']);
    %     print(figs(anNum).peak_diff,'-loose','-depsc2');
    
end






