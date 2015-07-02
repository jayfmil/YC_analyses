function YC1_subject_region_report(pool)


% analysis settings
% -----------------

% contrast to do
anas = {'correct_incorrect'};
ana_names = {'Correct - Incorrect'};

% limit to just correct trials?
only_correct = false;

% save directory
f = @(x,y) y{double(x)+1};
saveDirs = cell(1,length(anas));
for anNum = 1:length(anas)
    saveDirs{anNum} = fullfile('/Users/jmiller/Desktop/YC1/group',[anas{anNum},'_',f(only_correct,{'all_trials','correct_trials'})]);
    if ~exist(saveDirs{anNum},'dir')        
        error('Invalid analysis: %s',anas{anNum})
    elseif ~exist(fullfile(saveDirs{anNum},'figs'),'dir')
        mkdir(fullfile(saveDirs{anNum},'figs'));
    end
end

% get list of train subjects
% subjs = get_subs('RAM_YC1');
subjs = {'R1001P'}

if exist('pool','var')
    matlabpool(pool,length(subjs)+1)
    tic
    parfor s = 1:length(subjs)
        fprintf('Processing %s.\n',subjs{s})
        process_subj(subjs{s},anas,saveDirs,ana_names);
    end
    toc
    matlabpool close
    toc
else
    for s = subjs'
        fprintf('Processing %s.\n',s{1})
        process_subj(s{1},anas,saveDirs,ana_names);
    end
end

function process_subj(subj,anas,saveDirs,ana_names)

config = RAM_config('RAM_YC1');
rois = {'hipp'}%,'ec','mtl','frontal','parietal','temporal','occipital','limbic'};
figs = struct([]);

for anNum = 1:length(anas)
    close all

    
    saveDir = saveDirs{anNum};
    subj_files = dir(fullfile(saveDir,[subj,'*pow.mat']));
    if isempty(subj_files)
        fprintf('No files for analysis: %s %s\n', ana,subj)
        continue
    end
    ana = anas{anNum};
    figs(end+1).ana_name = ana_names{anNum};
    [~,ana_path] = fileparts(saveDir);
    figs(end).ana_path = ana_path;
     
    sig_countsLTA = NaN(length(rois),2);
    sig_countsHTA = NaN(length(rois),2);
    tstat_LTA = NaN(length(rois),2);
    tstat_HTA = NaN(length(rois),2);    

    sig_countsG = NaN(length(rois),2);
    sig_countsHFA = NaN(length(rois),2);
    tstat_G = NaN(length(rois),2);
    tstat_HFA = NaN(length(rois),2);


    sig_corrCountsLTA = NaN(length(rois),2);
    sig_corrCountsHTA = NaN(length(rois),2);
    freq_change = NaN(length(rois),1);
    freq_change_err = NaN(length(rois),1);
    num_elecs = NaN(length(rois),1);
    labels = {};
    link_text = {};
    
    figure('units','normalized','paperpositionmode','auto','position',[0.7504    0.6850    0.2188    0.2662])
    for r = 1:length(rois)
        
        if exist(fullfile(saveDir,[subj,'_',rois{r},'_pow.mat']),'file')
            regionData = load(fullfile(saveDir,[subj,'_',rois{r},'_pow.mat']));
        else
            labels{r} = '';
            continue
        end
        
        num_elecs(r) = length(regionData.statsLTA);
        labels{r} = sprintf('%s (%d)',rois{r},num_elecs(r));
        link_text{end+1} = rois{r};
        
        % low theta sig perc ttest
        ltaSig = [regionData.statsLTA.p] < .05;
        ltaCountPos = sum([regionData.statsLTA(ltaSig).tstat] > 0);
        sig_countsLTA(r,1) = ltaCountPos/num_elecs(r);
        sig_countsLTA(r,2) = (sum(ltaSig) - ltaCountPos)/num_elecs(r);        
        tstat_LTA(r,1) = nanmean([regionData.statsLTA.tstat]);
        tstat_LTA(r,2) = 1.96 * nanstd([regionData.statsLTA.tstat])/(sqrt(num_elecs(r))-1);        
        
        % high theta sig perc ttest
        htaSig = [regionData.statsHTA.p] < .05;
        htaCountPos = sum([regionData.statsHTA(htaSig).tstat] > 0);
        sig_countsHTA(r,1) = htaCountPos/num_elecs(r);
        sig_countsHTA(r,2) = (sum(htaSig) - htaCountPos)/num_elecs(r);        
        tstat_HTA(r,1) = nanmean([regionData.statsHTA.tstat]);
        tstat_HTA(r,2) = 1.96 * nanstd([regionData.statsHTA.tstat])/(sqrt(num_elecs(r))-1);         
        

        % gamma sig perc ttest
        gSig = [regionData.statsG.p] < .05;
        gCountPos = sum([regionData.statsG(gSig).tstat] > 0);
        sig_countsG(r,1) = gCountPos/num_elecs(r);
        sig_countsG(r,2) = (sum(gSig) - gCountPos)/ ...
            num_elecs(r);
        tstat_G(r,1) = nanmean([regionData.statsG.tstat]);
        tstat_G(r,2) = 1.96 * nanstd([regionData.statsG.tstat])/(sqrt(num_elecs(r))-1);

        % hfa sig perc ttest
        hfaSig = [regionData.statsHFA.p] < .05;
        hfaCountPos = sum([regionData.statsHFA(hfaSig).tstat] > 0);
        sig_countsHFA(r,1) = hfaCountPos/num_elecs(r);
        sig_countsHFA(r,2) = (sum(hfaSig) - hfaCountPos)/ ...
            num_elecs(r);
        tstat_HFA(r,1) = nanmean([regionData.statsHFA.tstat]);
        tstat_HFA(r,2) = 1.96 * ...
            nanstd([regionData.statsHFA.tstat])/(sqrt(num_elecs(r))-1);


        % low theta sig perc correlation
%        ltaSigCorr = [regionData.pLTA] < .05;
%        ltaCorrPos = sum([regionData.rLTA(ltaSigCorr)] > 0);
%        sig_corrCountsLTA(r,1) = ltaCorrPos/num_elecs(r);
%        sig_corrCountsLTA(r,2) = (sum(ltaSigCorr) - ltaCorrPos)/num_elecs(r);
        
        % high theta sig perc correlation
%        htaSigCorr = [regionData.pHTA] < .05;
%        htaCorrPos = sum([regionData.rHTA(htaSigCorr)] > 0);
%        sig_corrCountsHTA(r,1) = htaCorrPos/num_elecs(r);
%        sig_corrCountsHTA(r,2) = (sum(htaSigCorr) - htaCorrPos)/num_elecs(r);
        
        % frequency change
        freq_change(r) = nanmean(regionData.peakCond1 - regionData.peakCond2);
        freq_change_err(r) = nanstd(regionData.peakCond1 - regionData.peakCond2)/sqrt(num_elecs(r)-1) * 1.96;
        
        [h,p,c,s] = ttest([regionData.powCond1ByElec - regionData.powCond2ByElec]');
        powerDiff = nanmean(regionData.powCond1ByElec - regionData.powCond2ByElec,2)';
        powerDiff_STD = nanstd(regionData.powCond1ByElec - regionData.powCond2ByElec,[],2)/sqrt(num_elecs(r)-1);
        subplot(2,4,r)
        plot([log10(config.distributedParams.freQ(1)) log10(config.distributedParams.freQ(end))],[0 0],'-k','linewidth',1)
        hold on
        plotErrorRegions(log10(config.distributedParams.freQ),powerDiff,powerDiff_STD,{[.5 .5 .5]});
        hold on
        pthresh = .05/length(config.distributedParams.freQ);
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
        %     set(gca,'ylim',[-max(abs(get(gca,'ylim'))) max(abs(get(gca,'ylim')))])
        set(gca,'ylim',[-.8 .8])
%         set(gca,'xticklabel','')
        set(gca,'xticklabel',(2.^(0:ceil(sqrt(max(config.distributedParams.freQ))))))
        set(gca,'fontsize',8)
        title(labels{r},'fontsize',16)
        
        if r == 5
            set(gca,'xticklabel',(2.^(0:ceil(sqrt(max(config.distributedParams.freQ))))))
            xlabel('Frequency (Hz)','fontsize',14)
            ylabel('Z(Power) Diff.','fontsize',14)
        end
%         set(gca,'fontsize',14)
       
    end

    set(gcf,'PaperPositionMode','auto');
    figs(end).zpower = fullfile(saveDir,'figs',[subj '_' 'zpow_' ana '.eps']);    
    print(figs(end).zpower,'-depsc2');    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SIG ELECTRODE COUNT BAR PLOT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    clf
    colors = {[226,55,67]/255,[50,124,203]/255};
    h = bar(1:5:40,sig_countsLTA*100,.25,'stacked','linewidth',2);
    set(h,{'FaceColor'},{colors{1};colors{2}});
    hold on
    h = bar(2:5:41,sig_countsHTA*100,.25,'stacked','linewidth',2);
    set(h,{'FaceColor'},{colors{1};colors{2}});

    h = bar(3:5:42,sig_countsG*100,.25,'stacked','linewidth',2);
    set(h,{'FaceColor'},{colors{1};colors{2}});

    h = bar(4:5:43,sig_countsHFA*100,.25,'stacked','linewidth',2);
    set(h,{'FaceColor'},{colors{1};colors{2}});

    set(gca,'xtick',2.5:5:42)
    set(gca,'xticklabel',labels)
    set(gca,'ylim',[0 110]);
    set(gca,'ytick',[0:20:100]);
    rotateXLabels(gca,45)
    ylabel('% Sig. Elecs.','fontsize',16)
    set(gca,'fontsize',16)
    grid on
    set(gcf,'PaperPositionMode','auto');
    figs(end).sigCount_tTest = fullfile(saveDir,'figs',[subj '_' 'tTest_' ana '.eps']);    
    print(figs(end).sigCount_tTest,'-depsc2');
    keyboard
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SIG ELECTRODE CORR BAR PLOT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure(3)
%    clf
%    colors = {[226,55,67]/255,[50,124,203]/255};
%    h = bar(1:3:22,sig_corrCountsLTA*100,.25,'stacked','linewidth',2);
%    set(h,{'FaceColor'},{colors{1};colors{2}});
%    hold on
%    h = bar(2:3:23,sig_corrCountsHTA*100,.25,'stacked','linewidth',2);
%    set(h,{'FaceColor'},{colors{1};colors{2}});
%    set(gca,'xtick',1.5:3:23)
%    set(gca,'xticklabel',labels)
%    set(gca,'ylim',[0 110]);
%    set(gca,'ytick',[0:20:100]);
%    rotateXLabels(gca,45)
%    ylabel('% Sig. Elecs.','fontsize',16)
%    set(gca,'fontsize',16)
%    grid on
%    set(gcf,'PaperPositionMode','auto');
%    figs(end).sigCount_corr = fullfile(saveDir,'figs',[subj '_' 'corr_' ana '.eps']);    
%    print(figs(end).sigCount_corr,'-depsc2');    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% TSTAT BAR PLOT %%%
    %%%%%%%%%%%%%%%%%%%%%%

    figure(4)
    clf
    hold on
    colors = {[226,55,67]/255,[50,124,203]/255};    
    x1 = 1:5:40;
    x2 = 2:5:41;
    x3 = 3:5:42;
    x4 = 4:5:43;
    for i = 1:length(x1)
        c = colors{2};
        if tstat_LTA(i,1) > 0
            c = colors{1};
        end
        h=bar(x1(i),tstat_LTA(i,1),'linewidth',2);
        set(h,'FaceColor',c)
        errorbar(x1(i),tstat_LTA(i,1),tstat_LTA(i,2),'k','linewidth',2)

        c = colors{2};
        if tstat_HTA(i,1) > 0
            c = colors{1};
        end
        h=bar(x2(i),tstat_HTA(i,1),'linewidth',2);
        set(h,'FaceColor',c)
        errorbar(x2(i),tstat_HTA(i,1),tstat_HTA(i,2),'k','linewidth',2)        
      
        c = colors{2};
        if tstat_G(i,1) > 0
            c = colors{1};
        end
        h=bar(x3(i),tstat_G(i,1),'linewidth',2);
        set(h,'FaceColor',c)
        errorbar(x3(i),tstat_G(i,1),tstat_G(i,2),'k','linewidth',2)

  
        c = colors{2};
        if tstat_HFA(i,1) > 0
            c = colors{1};
        end
        h=bar(x4(i),tstat_HFA(i,1),'linewidth',2);
        set(h,'FaceColor',c)
        errorbar(x4(i),tstat_HFA(i,1),tstat_HFA(i,2),'k','linewidth',2)

    end

    set(gca,'xtick',2.5:5:42)
    set(gca,'xticklabel',labels)
    
    ylabel('Mean t-stat','fontsize',16)
    set(gca,'fontsize',16)
    grid on
    box on
%     if fixFiguresFor2014
%         set(gca,'XTickLabelRotation',-45)
%         set(gca,'GridLineStyle',':')
%         set(gca,'GridColor','k')
%         ylabel('Mean t-stat','fontsize',30)
%         set(gca,'fontsize',30)
%         set(gca,'clipping','off')
%     else
    rotateXLabels(gca,45)
%     end    
    
    set(gcf,'PaperPositionMode','auto');
    figs(end).tTest_tStat = fullfile(saveDir,'figs',[subj '_tTest_tstat_' ana '.eps']);    
    print(figs(end).tTest_tStat,'-depsc2');     
 keyboard   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PEAK FREQUENCY DIFF BAR PLOT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5)
    clf
    bar(1:length(rois),freq_change,'facecolor','w','linewidth',2);
    hold on
    errorbar(1:length(rois),freq_change,freq_change_err,'k','linestyle','none','linewidth',2)
    % set(gca,'xtick',1:1:length(rois)+1)
    set(gca,'xticklabel',labels)
    % set(gca,'ylim',[0 110]);
    % set(gca,'ytick',[0:20:100]);
    rotateXLabels(gca,45)
    set(gca,'xlim',[.25 length(rois)+.75]);
    ylabel('Peak Freq. Diff. (Hz)','fontsize',16)
    set(gca,'fontsize',16)
    grid on
    set(gcf,'PaperPositionMode','auto');
    figs(end).peak_diff = fullfile(saveDir,'figs',[subj '_' 'peak_diff_' ana '.eps']);    
    print(figs(end).peak_diff,'-depsc2');        
    keyboard
end


[baseDir,anaDir] = fileparts(saveDir);
texDir = fullfile(baseDir,'reports');
if ~exist(texDir,'dir')
    mkdir(texDir);
end

trials = '_all_trials';
if isempty(strfind(anaDir,'all_trials'))
    trials = '_correct_trials';
end
texName = [subj trials '_region_report.tex'];
write_texfile(texDir,texName,subj,figs,pthresh,link_text)


curr_dir = pwd;
cd(texDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(texDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);


% Start making the tex file
function write_texfile(saveDir,texName, subj, figs, pthresh, link_text)

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
fprintf(fid,'\\lhead{Train Subject Report: %s }\n',strrep(subj,'_','\_'));
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
for anNum = 1:length(figs)
    
%     fprintf(fid,'\\section*{%s}\n\n\n',figs(anNum).ana_name);    
    fprintf(fid,'\\href{run:%s}{\\Large\\bfseries %s}\n',fullfile('../..','pSpects',figs(anNum).ana_path,subj,[subj '_pSpect_report.pdf']),figs(anNum).ana_name);
    
    f = @(x) sprintf('\\href{%s#%s-%s}{\\small\\bfseries %s}',fullfile('../..','pSpects',figs(anNum).ana_path,subj,[subj '_pSpect_report.pdf']),subj,x,x);
    fprintf(fid,'\n\n\\small\\bfseries Links: ');
    for l = 1:length(link_text)
        fprintf(fid,'%s ',f(link_text{l}));
    end
    fprintf(fid,'\n');
    
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');       
    fprintf(fid,'\\includegraphics[width=0.99\\textwidth]{%s}\n',figs(anNum).zpower);
    fprintf(fid,'\\caption{\\textbf{Z score power difference by region.} Condition 1 - Condition 2. Positive values indicate greater power for condition 1. Significance threshold = %s.}\n',num2str(pthresh));
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\subfloat[]{\\includegraphics[width=0.43\\textwidth]{%s}}\n',figs(anNum).peak_diff);    
    fprintf(fid,'\\subfloat[]{\\includegraphics[width=0.43\\textwidth]{%s}}\n',figs(anNum).sigCount_tTest);   
    fprintf(fid,'\\subfloat[]{\\includegraphics[width=0.43\\textwidth]{%s}}\n',figs(anNum).tTest_tStat);
    
    if ~isempty(strfind(figs(anNum).ana_path,'fast'))
        fprintf(fid,'\\subfloat[]{\\includegraphics[width=0.43\\textwidth]{%s}}\n',figs(anNum).sigCount_corr);    

        fprintf(fid,[['\\caption{'] ...
            ['(a)  Difference in peak 1-10 Hz peak frequency by region for condition 1 peak - condition 2 peak, averaged across electrodes. '] ...
            ['(b)  Percent of electrodes with significant ($p<.05$) difference in power for between condition. Left: 1-4Hz. Right: 4-10 Hz. Red: positive effect. Blue: negative effect. ']...
            ['(c)  Mean t-stat, averaged across elecrodes, for LTA and HTA. Left: 1-4Hz. Right: 4-10 Hz.']...                        
            ['(d)  Percent of electrodes with significant ($p<.05$) correlation of power with movement speed. Left: 1-4Hz. Right: 4-10 Hz.']...
            ['}\n\n']]);
    else
        fprintf(fid,[['\\caption{'] ...
            ['(a)  Difference in peak 1-10 Hz peak frequency by region for condition 1 peak - condition 2 peak, averaged across electrodes. '] ...
            ['(b)  Percent of electrodes with significant ($p<.05$) difference in power for between condition. Left: 1-4Hz. Right: 4-10 Hz. Red: positive effect. Blue: negative effect. ']...            
            ['(c)  Mean t-stat, averaged across elecrodes, for LTA and HTA. Left: 1-4Hz. Right: 4-10 Hz.']...                        
            ['}\n\n']]);
    end
    fprintf(fid,'\\end{figure}\n\n\n');
    fprintf(fid,'\\clearpage\n\n\n');
end

fprintf(fid,'\\end{document}\n\n\n');










