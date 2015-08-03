function YC1_plotChanceDist(subjs,ana_dirs)


% analysis settings
% -----------------

if ~exist('ana_name','var') || isempty(ana_dirs)
    ana_dirs = {'CorrectedPower'};
    ana_dirs{2} = 'OrigPower';
    
    ana_names = {'Corrected Power'};
    ana_names{2} = 'Original Power';
    
end


saveDir = '/data10/scratch/jfm2/YC1/multi/reports_new';
if ~exist('saveDir','dir')
    mkdir(saveDir)
end

figDir = fullfile(saveDir,'figs');
if ~exist('saveDir','dir')
    mkdir(figDir)
end

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

figs = [];
perf_all = [];
perc_all = [];
quarts_all = [];
for s = 1:length(subjs)
    
    
    events = [];
    subj = subjs{s};
    figs_subj = [];
    figs_subj.subj = [];
    figs_subj.behavior = [];
    figs_subj.quarts = [];
    for a = 1:length(ana_dirs)
        
        figs_subj.(strrep(ana_names{a},' ','_')) = [];        
        
        dataDir = fullfile('/data10/scratch/jfm2/YC1/multi/lassoReg_allEncoding_binary_new',ana_dirs{a});
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
        
        
        chanceData = load(chanceFile);
        lassoData  = load(lassoFile);
        perc       = mean(lassoData.perf >= chanceData.perf_all);
        if a == 2
            perc_all = [perc_all;perc];
            perf_all = [perf_all;lassoData.perf];
            
            classProb = vertcat(lassoData.res.yPred{:});
            rec       = vertcat(lassoData.res.yTest{:});
            [classProbSort,ind] = sort(classProb);
            recSort = rec(ind);
            start = 1:(length(recSort)/4):length(recSort);
            stop = [start(2:end)-1 length(recSort)];
            
            quarts = [mean(recSort(start(1):stop(1))) mean(recSort(start(2):stop(2))) ...
                mean(recSort(start(3):stop(3))) mean(recSort(start(4):stop(4)))];
            quarts_all = [quarts_all;quarts];
            
            clf           
            h = bar(((quarts -.5)/.5)*100,'w','linewidth',2);            
            xlabel('Quartile of Classifier Estimate','fontsize',20)
            ylabel('Recall Change (%)','fontsize',20)
            set(gca,'fontsize',20)
            grid on
            if strcmp(subj,'R1059J_1')
                keyboard
            end
            fname = fullfile(figDir,[subj '_quart.eps']);
            figs_subj.quarts = fname;
            print('-depsc2','-loose',fname);
        end
        
        if isempty(events)
            clf
            events = get_sub_events('RAM_YC1',subj);
            events = events(strcmp({events.type},'NAV_TEST'));
            errs   = [events.respPerformanceFactor];
            [n,x] = hist(errs,25);
            bar(x,n/sum(n),1,'w','linewidth',2)
            xlabel('Performance Score','fontsize',20)
            ylabel('Prob.','fontsize',20)
            titleStr = sprintf('%s: median score = %.3f',subj,median(errs));
            title(strrep(titleStr,'_',' '),'fontsize',20)
            set(gca,'fontsize',20)
            set(gca,'xlim',[0 1]);
            fname = fullfile(figDir,[subj '_behaviorHist.eps']);
            figs_subj.behavior = fname;
            print('-depsc2','-loose',fname);
        end
        
        clf
        [n,x] = hist(chanceData.perf_all,25);
        bar(x,n/sum(n),1,'w','linewidth',2)
        xlim = get(gca,'xlim');
        
        set(gca,'xlim',[.5-max(abs(.5-xlim)) .5+max(abs(.5-xlim))])
        set(gca,'xlim',[.2 .8]);
        grid on
        xlabel('Classifier Performance','fontsize',20)
        ylabel('Prob.','fontsize',20)
        ylim = get(gca,'ylim');
        hold on
        plot([lassoData.perf lassoData.perf],[0 1],'--r','linewidth',3)
        set(gca,'ylim',ylim)
        titleStr = sprintf('%s %s: percentile = %.3f',subj,ana_names{a},perc);
        title(strrep(titleStr,'_',' '),'fontsize',20)
        set(gca,'fontsize',20)
        
        fname = fullfile(figDir,[subj '_' strrep(ana_names{a},' ','_') '.eps']);
        figs_subj.subj = subj;
        figs_subj.(strrep(ana_names{a},' ','_')) = fname;
        print('-depsc2','-loose',fname);
                                        
        
        
    end    
    
    events = [];
    figs = [figs;figs_subj];
end
keyboard

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
    fprintf(fid,'\\includegraphics[width=0.24\\textwidth]{%s}\n',figs(s).Original_Power);
    fprintf(fid,'\\includegraphics[width=0.24\\textwidth]{%s}\n',figs(s).Corrected_Power);
    fprintf(fid,'\\includegraphics[width=0.24\\textwidth]{%s}\n',figs(s).quarts);
    fprintf(fid,'\\includegraphics[width=0.24\\textwidth]{%s}\n',figs(s).behavior);
    fprintf(fid,'\\caption{%s}\n\n',strrep(figs(s).subj,'_',' '));
    fprintf(fid,'\\end{figure}\n\n\n');    
    if mod(s,3) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
end

fprintf(fid,'\\end{document}\n\n\n');






