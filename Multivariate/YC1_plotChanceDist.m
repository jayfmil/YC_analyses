function YC1_plotChanceDist(subjs,ana_dirs)


% analysis settings
% -----------------

if ~exist('ana_name','var') || isempty(ana_dirs)
    ana_dirs = {'CorrectedPower'};
    ana_dirs{2} = 'OrigPower';
    
    ana_names = {'Corrected Power'};
    ana_names{2} = 'Original Power';
    
end


saveDir = '/data10/scratch/jfm2/YC1/multi/reports';
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
sCount = 0;
for s = 1:length(subjs)
    
    subj = subjs{s};
    figs_subj = [];
    figs_subj.subj = [];
    for a = 1:length(ana_dirs)
        
        figs_subj.(strrep(ana_names{a},' ','_')) = [];
        
        dataDir = fullfile('/data10/scratch/jfm2/YC1/multi/lassoReg_allEncoding_binary',ana_dirs{a});
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
    figs = [figs;figs_subj];
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
    fprintf(fid,'\\includegraphics[width=0.45\\textwidth]{%s}\n',figs(s).Original_Power);
    fprintf(fid,'\\includegraphics[width=0.45\\textwidth]{%s}\n',figs(s).Corrected_Power);   
    fprintf(fid,'\\caption{%s}\n\n',strrep(figs(s).subj,'_',' '));
    fprintf(fid,'\\end{figure}\n\n\n');
    fprintf(fid,'\\end{figure}\n\n\n');
    if mod(s,2) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
end

fprintf(fid,'\\end{document}\n\n\n');






