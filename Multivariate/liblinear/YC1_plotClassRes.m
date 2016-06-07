function aucAll = YC1_plotClassRes(subjs,params,overwrite)
% function TH_plotClassRes(subjs,params,overwrite)
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
aucAll       = NaN(length(subjs),1);
recChangeAll = NaN(length(subjs),3);

% will hold figure paths for latex report
figs = [];
for s = 1:length(subjs)
    subj = subjs{s};
    fprintf('Creating plots for %s.\n',subj);
    
    % will hold subject specific figure paths
    figs_subj = struct('subj',[],'AUC',[],'recChange',[],'nElecs',[],'region',[]);
    if isempty(params.region)
        params.region = 'all';
    end
    figs_subj.subj   = subj;
    
    % see if files exist for subject. if not, continue
    classFile = fullfile(dataDir,[subj '_lasso_pow.mat']);
    if params.usePhase==1
        classFile = fullfile(dataDir,[subj '_lasso_phase.mat']);
    elseif params.usePhase==2
        classFile = fullfile(dataDir,[subj '_lasso_powphase.mat']);
    end
    if ~exist(classFile,'file')
        fprintf('Lasso not found for %s.\n',subj)
        continue
    end
    
    % load subject data
    classData  = load(classFile);
    
    % store informatiopn about subject
    figs_subj.nElecs   = length(classData.tal);
    figs_subj.region   = classData.params.region;    
    if isempty(figs_subj.region);figs_subj.region='All';end
    nTrialsAll(s)      = length(classData.Y);
    events             = get_sub_events('RAM_YC1',subj);
    
    % store info for group plots
    aucAll(s,:)    = classData.AUC;
    
    %--------------------------------------------------------------------------
    % Figure 1 - recall chance as a function of classifier output plot
    
    % sort recalled or not by classifier output
    rec       = vertcat(classData.res.yTest{:});
    classProb = vertcat(classData.res.yPred{:});
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
        figure(1)
        clf
        h = bar(recChangeBin,'w','linewidth',2);
        xlabel('Tercile of Classifier Estimate','fontsize',20)
        ylabel('Recall Change (%)','fontsize',20)
        set(gca,'fontsize',20)
        y = get(gca,'ylim');
        set(gca,'ylim',[-1 1]*max(abs(y)))
        set(gca,'xlim',[0 4])
        grid on
        set(gca,'gridlinestyle',':');
        hold on
        print('-depsc2','-tiff','-loose',fname);
    end
    
    
    %--------------------------------------------------------------------------
    % Figure 2 - ROC
    
    fname = fullfile(figDir,[subj '_AUC.eps']);
    figs_subj.AUC = fname;
    if ~exist(fname,'file') || overwrite
        
        figure(2)
        clf
        labels = vertcat(classData.res.yTest{:});
        scores = vertcat(classData.res.yProb{:});
        [x,y] = perfcurve(labels,scores,1);
        plot(x,y,'-k','linewidth',3)
        hold on
        plot([0 1],[0 1],'-k','linewidth',2,'color',[.5 .5 .5])
        grid on
        set(gca,'gridlinestyle',':')
        xlabel('False Alarm Rate','fontsize',20);
        ylabel('Hit Rate','fontsize',20);
        axis square
        set(gca,'xtick',0:.2:1)
        set(gca,'ytick',0:.2:1)
        set(gca,'fontsize',20);
        titleStr = sprintf('AUC = %.3f',classData.res.AUC);
        title(titleStr)
        print('-depsc2','-tiff','-loose',fname);
    end
    figs = [figs;figs_subj];        
end

good = ~cellfun('isempty',{figs.subj});
figs = figs(good);
texName = 'classReport.tex';
write_texfile(saveDir,texName,figs)


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
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).recChange);
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).AUC);
    fprintf(fid,'\\caption{%s: region: %s, %d electrodes}\n\n',strrep(figs(s).subj,'_',' '),figs(s).region,figs(s).nElecs);
    fprintf(fid,'\\end{figure}\n\n\n');
    if mod(s,2) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
end

fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);









