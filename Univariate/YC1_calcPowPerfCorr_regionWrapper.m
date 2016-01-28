function YC1_calcPowPerfCorr_regionWrapper(subjs)


% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% regions to iterate over
regions = {'Hippform','Hipp','CA1','CA2','CA3','DG','Sub','EC'};

params = univarParams();
params.basePath = '/scratch/jfm2/YC1/uniPreAveragedPower/learnTrials';
params.overwrite = 0;
% for each region, do both averaged and non-averaged over electrodes
figs = [];
figsPower = [];
for r = 1:length(regions)
    figs(1).(regions{r}) = struct();
    
    params.region = regions{r};
    params.averageRegion = 0;
    YC1_calcPowPerfCorr_usePreAveragedPower(subjs,params)
    figs(1).(regions{r}).allElecs = YC1_makeCorrPlots(subjs,params);
    figsPower(1).(regions{r}).allElecs = YC1_plotPowerSpect(subjs,params); 
    params.averageRegion = 1;
    YC1_calcPowPerfCorr_usePreAveragedPower(subjs,params);
    figs(1).(regions{r}).avgElecs = YC1_makeCorrPlots(subjs,params);
    figsPower(1).(regions{r}).avgElecs = YC1_plotPowerSpect(subjs,params); 
end

saveDir = fullfile(params.basePath,'OrigPower','report');
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

texName = 'univariateLearningReport.tex';
write_texfile(saveDir,texName,figs)

curr_dir = pwd;
cd(saveDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(saveDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);

texName = 'powSpectLearningReport.tex';
write_powSpect_texfile(saveDir,texName,figs)

curr_dir = pwd;
cd(saveDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(saveDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);


params.basePath = '/scratch/jfm2/YC1/uniPreAveragedPower/testTrials';
params.eventFilter=@(events)(strcmp({events.type},'NAV_TEST')) ;
% for each region, do both averaged and non-averaged over electrodes
figs = [];
for r = 1:length(regions)
    figs(1).(regions{r}) = struct();
    
    params.region = regions{r};
    params.averageRegion = 0;
    YC1_calcPowPerfCorr_usePreAveragedPower(subjs,params)
    figs(1).(regions{r}).allElecs = YC1_makeCorrPlots(subjs,params);
    params.averageRegion = 1;
    YC1_calcPowPerfCorr_usePreAveragedPower(subjs,params);
    figs(1).(regions{r}).avgElecs = YC1_makeCorrPlots(subjs,params);
    
end

saveDir = fullfile(params.basePath,'OrigPower','report');
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end


texName = 'univariateTestReport.tex';
write_texfile(saveDir,texName,figs)

curr_dir = pwd;
cd(saveDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(saveDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);

texName = 'powSpectTestReport.tex';
write_powSpect_texfile(saveDir,texName,figs)

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
fprintf(fid,'\\usepackage{subfigure,amsmath} \n');
% fprintf(fid,'\\usepackage{wrapfig} \n');
% fprintf(fid,'\\usepackage{longtable} \n');
% fprintf(fid,'\\usepackage{pdfpages}\n');
% fprintf(fid,'\\usepackage{mathtools}\n');
% fprintf(fid,'\\usepackage{array}\n');
% fprintf(fid,'\\usepackage{enumitem}\n');
% fprintf(fid,'\\usepackage{sidecap} \\usepackage{soul}\n');

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
fprintf(fid,'\\lhead{All Subjects}\n');
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
fields = fieldnames(figs);

for f = 1:length(fields)
    
    allElecs = figs.(fields{f}).allElecs;
    fprintf(fid,'\\begin{figure}\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.33\\textwidth]{%s}}}\n',allElecs.corrByFreq);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.33\\textwidth]{%s}}}\n',allElecs.corrByFreqInner);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.33\\textwidth]{%s}}}\n',allElecs.corrByFreqOuter);    
    fprintf(fid,'\\caption{Region %s, %d electrodes average. a: All locations. b: Inner rectangle. c: Outer region.}\n',fields{f},allElecs.numElecs);
    fprintf(fid,'\\end{figure}\n\n\n');
    
    avgElecs = figs.(fields{f}).avgElecs;
    fprintf(fid,'\\begin{figure}\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.33\\textwidth]{%s}}}\n',avgElecs.corrByFreq);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.33\\textwidth]{%s}}}\n',avgElecs.corrByFreqInner);
    fprintf(fid,'\\subfigure[]{{\\includegraphics[width=0.33\\textwidth]{%s}}}\n',avgElecs.corrByFreqOuter);    
    fprintf(fid,'\\caption{Region %s, %d subjects average. a: All locations. b: Inner rectangle. c: Outer region.}\n',fields{f},avgElecs.numElecs);
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\clearpage\n\n\n');
    
end

fprintf(fid,'\\end{document}\n\n\n');




% Start making the tex file
function write_powSpect_texfile(saveDir,texName, figs)

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
fprintf(fid,'\\usepackage{subfigure,amsmath} \n');
fprintf(fid,'\\usepackage{longtable} \n');
% fprintf(fid,'\\usepackage{wrapfig} \n');
% fprintf(fid,'\\usepackage{longtable} \n');
% fprintf(fid,'\\usepackage{pdfpages}\n');
% fprintf(fid,'\\usepackage{mathtools}\n');
% fprintf(fid,'\\usepackage{array}\n');
% fprintf(fid,'\\usepackage{enumitem}\n');
% fprintf(fid,'\\usepackage{sidecap} \\usepackage{soul}\n');

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
fprintf(fid,'\\lhead{All Subjects}\n');
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
fields = fieldnames(figs);

for f = 1:length(fields)
        
    allElecs = figs.(fields{f}).allElecs;
    subjs = fieldnames(allElecs);
    for s = 1:length(subjs)
        subj = subjs{s};
        fprintf(fid,'\\section{%s: %s}\n',strrep(subj,'_',' '),fields{f});
        fprintf(fid,'\\begin{longtable}{cc}\n');
        
        subjFigs = allElecs.(subj);
        for elecNum = 1:length(subjFigs)                        
            fprintf(fid,'\\includegraphics[width=.45\\textwidth]{%s}',subjFigs(elecNum).fname);
            if mod(elecNum,2) == 1
                fprintf(fid,'&');
            else
                fprintf(fid,'\\\\\n');
            end
        end
        if elecNum >= 1 && mod(elecNum,2) == 1
            fprintf(fid,'\\\\\n');
        end        
        fprintf(fid,'\\end{longtable}\n');        
    end
    
    
end

fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);





