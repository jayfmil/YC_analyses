function YC1_makeUnivarSubjectReports(subjs,params,overwrite)
% function YC1_makeUnivarSubjectReports(subjs,params,overwrite)
%
% Inputs:
%
%      subjs - cell array of subject strings (default get_subs('RAM_YC1'))
%     params - params structure (default is returned by univarParams)
%  overwrite - boolean. if true, overwrite existing figures
%
%
%
%
%

% if not given, use default params
if ~exist('params','var') || isempty(params)
    params = univarParams();
end

% overwrite exisiting figures?
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

% save directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
region = params.region;
if isempty(region);region = 'all';end
dataDir = fullfile(params.basePath,f(params.useCorrectedPower,y),region);
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
nTrialsAll = NaN(length(subjs),1);
perf_all   = NaN(length(subjs),size(params.timeBins,1));
perf_p_all = NaN(length(subjs),size(params.timeBins,1));
auc_all    = NaN(length(subjs),size(params.timeBins,1));
auc_p_all  = NaN(length(subjs),size(params.timeBins,1));
quarts_all = NaN(size(params.timeBins,1),4,length(subjs));
best_time  = NaN(length(subjs),1);

% will hold figure paths for latex report
figs = [];
for s = 1:length(subjs)                
    subj = subjs{s};
    fprintf('Creating plots for %s.\n',subj);
        
    % will subject specific figure paths
    figs_subj = struct('subj',[],'behavior',[],'quarts',[],'region','',...
                       'AUC','','perf','','nElecs',[]);        
    if isempty(params.region)
        params.region = 'all';
    end
    figs_subj.region = params.region;
    figs_subj.subj   = subj;
    
    % see if files exist for subject. if not, continue
    chanceFile = fullfile(dataDir,[subj '_chance_perf_dist.mat']);
    lassoFile  = fullfile(dataDir,[subj '_lasso.mat']);    
    if ~exist(chanceFile,'file')
        fprintf('Chance distribution not found for %s.\n',subj)
        continue
    end
    if ~exist(lassoFile,'file')
        fprintf('Lasso not found for %s.\n',subj)
        continue
    end
    
    % load lasso and chance lasso data
    chanceData = load(chanceFile);
    lassoData  = load(lassoFile);
    if length(lassoData.Y) < 48
        fprintf('Less than half a session of data for %s. Skipping.\n',subj)
        continue
    end
    
    figs_subj.nElecs = length(lassoData.res(1).A{1})/4;
    nTrialsAll(s) = length(lassoData.Y);            
    events     = get_sub_events('RAM_YC1',subj);    
    

    
    % calculate percentile for subject
    perf_p          = mean(repmat(lassoData.perf,size(chanceData.perf_all,1),1) > chanceData.perf_all);
    auc_p           = mean(repmat(lassoData.AUC,size(chanceData.auc_all,1),1) > chanceData.auc_all);
    perf_p_all(s,:) = perf_p;
    auc_p_all(s,:)  = auc_p;
    perf_all(s,:)   = lassoData.perf;
    auc_all(s,:)    = lassoData.AUC;
    
    
    %% FIGURE 1a and b - classifier accuracy and AUC over time       
    
    % create x labels based on time bins
    xBins    = round((params.timeBins) / 10) / 100;
    xBinsStr = {};
    for x = 1:size(xBins,1)
        xBinsStr{x} = [num2str(xBins(x,1)), '-', num2str(xBins(x,2))];
    end
        
    ylabels = {'Classifier Accuracy (%)','Classifier AUC'};
    ps      = [1-perf_p;1-auc_p];
    fields  = {'perf','AUC'};

    f = @(x) nnz(mean(horzcat(lassoData.res(x).A{:}),2));
    emptyModels = cellfun(f,num2cell(1:length(lassoData.res)))==0;

    for i = 1:2
        fname = fullfile(figDir,[subj '_' fields{i} '.eps']);
        figs_subj.(fields{i}) = fname;
        if exist(fname,'file') && ~overwrite
            continue
        end
        figure(1) 
        clf
        
        % first plot all the points as black
        plot(lassoData.(fields{i}),'k.','markersize',30,'linewidth',2.5)
        hold on

        % plot all empty models as X
        plot(find(emptyModels),lassoData.(fields{i})(emptyModels),'xk','markersize',65)
        
        % plot any significant time points as red. Bigger red dot indicates
        % the timepoint survived bonferroni correction. Small red dot is p
        % < .05        
        p = ps(i,:);
        thresh = .05/length(lassoData.(fields{i}));
        h=plot(find(p < .05),lassoData.(fields{i})(p < .05),'r.','markersize',30);
        set(h,'color',[140 15 15]/255)
        plot(find(p < thresh),lassoData.(fields{i})(p < thresh),'r.','markersize',55)
        set(h,'color',[200 100 100]/255)
        
        % set axis and labels
        set(gca,'xtick',1:length(lassoData.(fields{i})));
        set(gca,'xlim',[0 length(lassoData.(fields{i}))+1]);
        set(gca,'xticklabel',xBinsStr);
        set(gca,'ylim',[0 1]);
        set(gca,'ytick',0:.1:1);
        if i == 1
            set(gca,'yticklabel',0:10:100)
        end
        ylabel(ylabels{i},'fontsize',20)
        xlabel('Time (s)','fontsize',20)
        grid on
        set(gca,'gridlinestyle',':');
        set(gca,'fontsize',20)
        fillX = [.5 .5 length(lassoData.(fields{i}))+.05 length(lassoData.(fields{i}))+.05];
        h=fill(fillX,[.91 .98 .98 .91],'w');
        set(h,'linestyle','none');
        
        % label each time point with the behavior
        labels = params.timeBinLabels;
        if isempty(labels);labels=repmat({''},1,size(params.timeBins,1));end
        for t = 1:length(lassoData.(fields{i}))
            start = t - .25;
            stop  = t + .25;
            plot([start stop],[.9 .9],'-k','linewidth',2)
            h=text(t,.95,labels{t},'fontsize',16);
            len = get(h,'extent');
            len = len(3);
            pos = get(h,'position');
            set(h,'position',[pos(1)-len/2+.03 pos(2) 0]);
        end
        
        % plot dashed line at 50%
        xlim = get(gca,'xlim');
        plot(xlim,[.5 .5],'--k','linewidth',1.5)
        box on        
                        
        print('-depsc2','-tiff','-loose',fname);
    end
    
    %% FIGURE 2 - quartile plot for most signficant time period
    fname = fullfile(figDir,[subj '_quart.eps']);
    figs_subj.quarts = fname;
    
    % recalled vs not recalled vector
    rec       = vertcat(lassoData.res(1).yTest{:});
    
    % compute quartiles
    bestTime  = find(perf_p == max(perf_p),1,'first');
    best_time(s) = bestTime;
    for t = 1:length(lassoData.res)
        classProb = vertcat(lassoData.res(t).yPred{:});
        
        % sort the recall vector based on the classifier probs
        [~,ind] = sort(classProb);
        recSort = rec(ind);
        
        % now bin the sorted recall vector
        start = 1:(length(recSort)/4):length(recSort);
        stop = [start(2:end)-1 length(recSort)];
        quarts_all(t,:,s) = [mean(recSort(start(1):stop(1))) mean(recSort(start(2):stop(2))) ...
            mean(recSort(start(3):stop(3))) mean(recSort(start(4):stop(4)))];
        if t == bestTime
            quarts = quarts_all(t,:,s);
        end
    end
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(2)
        clf
        % plot quartiles based on change from mean
        h = bar(((quarts -(mean(rec)))/mean(rec))*100,'w','linewidth',2);
        xlabel('Quartile of Classifier Estimate','fontsize',20)
        ylabel('Recall Change (%)','fontsize',20)
        set(gca,'fontsize',20)
        set(gca,'ylim',[-100 100])
        set(gca,'xlim',[0 5])
        grid on
        set(gca,'gridlinestyle',':');
        hold on
        h=title([labels{bestTime} ' Period'],'fontsize',20);
        set(h,'fontweight','normal');
        
        print('-depsc2','-tiff','-loose',fname);
    end
    
    %% Figure 3 - response error histogram
    fname = fullfile(figDir,[subj '_behavior.eps']);
    figs_subj.behavior = fname;
    if exist(fname,'file') && overwrite || (~exist(fname,'file'))
        figure(3)
        clf
        
        % getting this from the events structure. Really i should save it in
        % the results of the classifier as an unbinarized version of Y, in case
        % we use a difference performance measure
        errs = [events(strcmp({events.type},'NAV_TEST')).respPerformanceFactor];
        
        % plot histogram
        [n,x] = hist(errs,25);
        bar(x,n/sum(n),1,'w','linewidth',2)
        xlabel('Performance Score','fontsize',20)
        ylabel('Prob.','fontsize',20)
        titleStr = sprintf('%s: median score = %.3f',subj,median(errs));
        h=title(strrep(titleStr,'_',' '),'fontsize',20);
        set(h,'fontweight','normal') ;
        set(gca,'fontsize',20)
        set(gca,'xlim',[0 1]);
        
        print('-depsc2','-tiff','-loose',fname);
    end
    figs = [figs;figs_subj];
    
end

% also make group plots/report
fprintf('Creating group plots.\n');
figs_group = [];
figs_group.quarts   = {};
figs_group.acc_hist = {};
figs_group.auc_hist = {};

% compute average quartile measure for each time bin
meanRec_subj = repmat(nanmean(quarts_all,2),1,4,1);
nSubj        = sum(~isnan(quarts_all),3);
figs_group.N = nSubj(:,1);
quarts_err   = nanstd((quarts_all - meanRec_subj)./meanRec_subj,[],3)./sqrt(nSubj-1);
quarts_group = nanmean((quarts_all - meanRec_subj)./meanRec_subj,3);

% labels for plotting
labels = params.timeBinLabels;
if isempty(labels);labels=repmat({''},1,size(quarts_group,1));end

% plot each time bin separately
for t = 1:size(quarts_group,1)
    
    %% QUARTILE PLOT
    fname = fullfile(figDir,['group_quart_' labels{t} '.eps']);
    figs_group.quarts{t} = fname;
    
    if (exist(fname,'file') && overwrite) || (~exist(fname,'file'))
        figure(2)
        clf
        bar(quarts_group(t,:)*100,'w','linewidth',2);
        hold on
        errorbar(1:4,quarts_group(t,:)*100,quarts_err(t,:)*196,'k','linewidth',2,'linestyle','none')
        
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
fprintf(fid,'\\caption{%d Subjects: Subject average quartile by time bin.}\n\n',figs.N(1,1));
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



function sout = mergestruct(varargin)
%MERGESTRUCT Merge structures with unique fields.
%   Copyright 2009 The MathWorks, Inc.

% Start with collecting fieldnames, checking implicitly
% that inputs are structures.
fn = [];
for k = 1:nargin
    try
        fn = [fn ; fieldnames(varargin{k})];
    catch MEstruct
        throw(MEstruct)
    end
end

% Make sure the field names are unique.
if length(fn) ~= length(unique(fn))
    error('mergestruct:FieldsNotUnique',...
        'Field names must be unique');
end

% Now concatenate the data from each struct.  Can't use
% structfun since input structs may not be scalar.
c = [];
for k = 1:nargin
    try
        c = [c ; struct2cell(varargin{k})];
    catch MEdata
        throw(MEdata);
    end
end
% Construct the output.
sout = cell2struct(c, fn, 1);



