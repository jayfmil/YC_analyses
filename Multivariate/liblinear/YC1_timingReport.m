function YC1_timingReport(subjs,params,overwrite)
% function YC1_timingReport(subjs,params,overwrite)
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
%   - Classifier AUC over time
%   - Recall change as a function of classifier output for the best
%     performing  time bin
%
%
% For the group report, figures are Tercile plots for each time bin,
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
dataDir = fullfile(params.basePath,f(params.useCorrectedPower,y));
saveDir = fullfile(dataDir,'reports');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% figure directory
figDir = fullfile(saveDir,'figs');
if ~exist('figDir','dir')
    mkdir(figDir)
end

% store group information
nTrialsAll   = NaN(length(subjs),1);
aucAll       = NaN(length(subjs),size(params.timeBins,1));
recChangeAll = NaN(length(subjs),3);
bestTimeAll  = NaN(length(subjs),1);

% will hold figure paths for latex report
figs = [];
for s = 1:length(subjs)                
    subj = subjs{s};
    fprintf('Creating plots for %s.\n',subj);
        
    % will hold subject specific figure paths
    figs_subj = struct('subj',[],'AUC',[],'recChange',[],'nElecs',[],'bestTime',[]);        
    if isempty(params.region)
        params.region = 'all';
    end
    figs_subj.subj   = subj;
    
    % see if files exist for subject. if not, continue    
    lassoFile  = fullfile(dataDir,[subj '_lasso.mat']);    
    if ~exist(lassoFile,'file')
        fprintf('Lasso not found for %s.\n',subj)
        continue
    end
    
    % load subject data    
    lassoData  = load(lassoFile);    
    
    % store informatiopn about subject
    bestTime           = find(lassoData.AUC == max(lassoData.AUC),1,'first');
    figs_subj.bestTime = bestTime;
    figs_subj.nElecs   = length(lassoData.tal);
    nTrialsAll(s)      = length(lassoData.Y);            
    events             = get_sub_events('RAM_YC1',subj);    
    
    % store info for group plots
    aucAll(s,:)    = lassoData.AUC;
    bestTimeAll(s) = bestTime;
%     bestTimeAll(s) = 9;
    
%--------------------------------------------------------------------------
% Figure 1 - classifier AUC over time
       
    fname = fullfile(figDir,[subj '_AUC_byTime.eps']);
    figs_subj.AUC = fname;
    if ~exist(fname,'file') || overwrite
        
        % create time bin labels
        xBins    = round((params.timeBins) / 10) / 100;
        xBinsStr = {};
        for x = 1:size(xBins,1)
            xBinsStr{x} = [num2str(xBins(x,1)), '-', num2str(xBins(x,2))];
        end
        
        % see if any times bins have empty models, will be marked by X on
        % plot
        f = @(x) nnz(mean(horzcat(lassoData.res(x).A{:}),2));
        emptyModels = cellfun(f,num2cell(1:length(lassoData.res)))==0;
        
        % plot it
        figure(1)
        clf
        plot(lassoData.AUC,'-k','linewidth',3);
        ylabel('AUC','fontsize',20)
        xlabel('Time (s)','fontsize',20)
        grid on
        set(gca,'xlim',[0 length(lassoData.AUC)+1])
        set(gca,'xtick',1:10:length(lassoData.AUC));
        set(gca,'xticklabel',xBinsStr(1:10:length(lassoData.AUC)))
        set(gca,'gridlinestyle',':');
        set(gca,'fontsize',20)
        set(gca,'xticklabelrotation',45) 
        hold on
        plot(find(emptyModels),lassoData.AUC(emptyModels),'xk','markersize',65)
        plot([0 length(lassoData.AUC)+1],[.5 .5],'--k','linewidth',2)
        ylim = get(gca,'ylim');
        plot([bestTime bestTime],ylim,'--r','linewidth',2)
        set(gca,'ylim',ylim)
        print('-depsc2','-tiff','-loose',fname);         
    end
    
%--------------------------------------------------------------------------    
% Figure 2 - recall chance as a function of classifier output plot

    % sort recalled or not by classifier output
    rec       = vertcat(lassoData.res(bestTime).yTest{:});
    classProb = vertcat(lassoData.res(bestTime).yPred{:});
    [~,ind] = sort(classProb);
    recSort = rec(ind);
    recSort(recSort==-1) = 0;
    meanRec = mean(recSort);
    
    % bin into terciles
    bins  = [0 round(length(recSort)/3) * [1:2] length(recSort)];
    start = bins(1:end-1) + 1;
    stop  = bins(2:end);
    recChangeBin = NaN(1,length(start));
    for i = 1:length(start)
        recChangeBin(i) = mean(recSort(start(i):stop(i)));
    end
    recChangeBin = (recChangeBin-meanRec)/meanRec * 100;
    
    % store for group plot
    recChangeAll(s,:) = recChangeBin;
    
    fname = fullfile(figDir,[subj '_recChange.eps']);
    figs_subj.recChange = fname;    
    if ~exist(fname,'file') || overwrite
        
        % plot it
        figure(2)
        clf        
        h = bar(recChangeBin,'w','linewidth',2);
        xlabel('Tercile of Classifier Estimate','fontsize',20)
        ylabel('Recall Change (%)','fontsize',20)
        set(gca,'fontsize',20)
        set(gca,'ylim',[-40 40])
        set(gca,'xlim',[0 4])
        grid on
        set(gca,'gridlinestyle',':');
        hold on       
        print('-depsc2','-tiff','-loose',fname);        
    end        
    figs = [figs;figs_subj];
    
end


%--------------------------------------------------------------------------    
% Figure 3 - group AUC over time
figs_group = [];
mAuc = nanmean(aucAll);
sAuc = nanstd(aucAll);
eAUC = sAuc./sqrt(sum(~isnan(aucAll))-1);

fname = fullfile(figDir,['group_AUC_byTime.eps']);
figs_group.AUC = fname;
if ~exist(fname,'file') || overwrite
    
    % create time bin labels
    xBins    = round((params.timeBins) / 10) / 100;
    xBinsStr = {};
    for x = 1:size(xBins,1)
        xBinsStr{x} = [num2str(xBins(x,1)), '-', num2str(xBins(x,2))];
    end
    
    % plot it
    figure(3)
    clf
    plot(mAuc,'-k','linewidth',3);
    ylabel('AUC','fontsize',20)
    xlabel('Time (s)','fontsize',20)
    grid on
    set(gca,'xlim',[0 length(mAuc)+1])
    set(gca,'xtick',1:10:length(mAuc));
    set(gca,'xticklabel',xBinsStr(1:10:length(mAuc)))
    set(gca,'gridlinestyle',':');
    set(gca,'fontsize',20)
    set(gca,'xticklabelrotation',45)
    hold on    
    plot([0 length(mAuc)+1],[.5 .5],'--k','linewidth',2)
    ylim = get(gca,'ylim');
%     plot([bestTime bestTime],ylim,'--r','linewidth',2)
    set(gca,'ylim',ylim)
    print('-depsc2','-tiff','-loose',fname);
end


%--------------------------------------------------------------------------    
% Figure4 - group tercile
mRec = nanmean(recChangeAll);
sRec = nanstd(recChangeAll);
eRec = sRec./sqrt(sum(~isnan(recChangeAll))-1);

fname = fullfile(figDir,['group_recChange.eps']);
figs_group.recChange = fname;
if ~exist(fname,'file') || overwrite
    
    % plot it
    figure(4)
    clf
    bar(mRec,'w','linewidth',2);
    hold on
    errorbar(1:3,mRec,eRec*1.96,'-k','linewidth',2,'linestyle','none')
    xlabel('Tercile of Classifier Estimate','fontsize',20)
    ylabel('Recall Change (%)','fontsize',20)
    set(gca,'fontsize',20)
    set(gca,'ylim',[-100 100])
    set(gca,'xlim',[0 4])
    grid on
    set(gca,'gridlinestyle',':');
    
    print('-depsc2','-tiff','-loose',fname);
end
keyboard

% inds = sub2ind(size(aucAll),[1:length(bestTimeAll)]',bestTimeAll)

% also make group plots/report
fprintf('Creating group plots.\n');
figs_group = [];
figs_group.quarts   = {};
figs_group.acc_hist = {};
figs_group.auc_hist = {};

% compute average Tercile measure for each time bin
meanRec_subj = repmat(nanmean(quarts_all,2),1,3,1);
nSubj        = sum(~isnan(quarts_all),3);
figs_group.N = nSubj(:,1);
quarts_err   = nanstd((quarts_all - meanRec_subj)./meanRec_subj,[],3)./sqrt(nSubj-1);
quarts_group = nanmean((quarts_all - meanRec_subj)./meanRec_subj,3);

% labels for plotting
labels = params.timeBinLabels;
if isempty(labels);labels=repmat({''},1,size(quarts_group,1));end

% plot each time bin separately
for t = 1:size(quarts_group,1)
    
    %% Tercile PLOT
    fname = fullfile(figDir,['group_quart_' labels{t} '.eps']);
    figs_group.quarts{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(2)
        clf
        bar(quarts_group(t,:)*100,'w','linewidth',2);
        hold on
        errorbar(1:3,quarts_group(t,:)*100,quarts_err(t,:)*196,'k','linewidth',2,'linestyle','none')
        
        xlabel('Tercile of Classifier Estimate','fontsize',20)
        ylabel('Recall Change (%)','fontsize',20)
        set(gca,'fontsize',20)
        set(gca,'ylim',[-100 100])
        set(gca,'xlim',[0 4])
        set(gca,'ylim',[-25 25])
        grid on
        set(gca,'gridlinestyle',':');
        hold on
        h=title([labels{t} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');
        print('-depsc2','-tiff','-loose',fname);
    end
    
    %% ACCURACY HISTOGRAM PLOT
    fname = fullfile(figDir,['acc_hist_' labels{t} '.eps']);
    figs_group.acc_hist{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(3)
        clf
        perf = perf_all(:,t);
        sig  = perf_p_all(:,t) > .95;
        n1   = histc(perf(sig),0.025:.05:.975);
        % wtf
        if isrow(n1);n1=n1';end        
        n2   = histc(perf(~sig),0.025:.05:.975);
        h    = bar([.05:.05:1]*100,[n1 n2],1,'stacked','linewidth',2);
        xlabel('Classifier Percent Correct','Fontsize',20);
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
    
    %% AUC HISTOGRAM PLOT
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


good = ~cellfun('isempty',{figs.subj});
figs = figs(good);
texName = 'lassoChance_report.tex';
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
texName = 'group_lassoChance_report.tex';
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
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).perf);
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).AUC);
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).quarts);
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).behavior);    
    fprintf(fid,'\\caption{%s: region: %s, %d electrodes}\n\n',strrep(figs(s).subj,'_',' '),figs(s).region,figs(s).nElecs);
    fprintf(fid,'\\end{figure}\n\n\n');
    if mod(s,2) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
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
fprintf(fid,'\\caption{%d Subjects: Subject average Tercile by time bin.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');

fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)       
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.acc_hist{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject accuracy histogram by time bin.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');

fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
for f = 1:size(figs.N,1)       
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs.auc_hist{f});
end
fprintf(fid,'\\caption{%d Subjects: Subject AUC histogram by time bin}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');

fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');

fprintf(fid,'\\includegraphics[width=0.8\\textwidth]{%s}\n',figs.auc_time_direct);
fprintf(fid,'\\caption{%d Subjects: Average AUC over time. Dark red: significant after correcting for number of time bins. Light red: $p<.05$.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');

fprintf(fid,'\\end{document}\n\n\n');







