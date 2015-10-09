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
% nTrialsAll = NaN(length(subjs),1);
% perf_all   = NaN(length(subjs),size(params.timeBins,1));
% perf_p_all = NaN(length(subjs),size(params.timeBins,1));
% auc_all    = NaN(length(subjs),size(params.timeBins,1));
% auc_p_all  = NaN(length(subjs),size(params.timeBins,1));
% quarts_all = NaN(size(params.timeBins,1),4,length(subjs));
% best_time  = NaN(length(subjs),1);

% will hold figure paths for latex report
figs = [];
subjDataAll =  [];
for s = 1:length(subjs)                
    
    subj = subjs{s};
    subjFile = fullfile(dataDir,[subj '.mat']);
    if ~exist(subjFile,'file')
        fprintf('Subject data not found for %s.\n',subj)
        continue
    end    
    

    % load subejct data
    subjData = load(subjFile);
    
    % if we haven't done it yet, initialze structure to hold concatenated
    % data from all subjects
    if isempty(subjDataAll)
       fields = fieldnames(subjData.res);
       for f = fields'           
          subjDataAll.(f{1}) = []; 
          subfields = fieldnames(subjData.res.(f{1}));
          for subf = subfields'
              subjDataAll.(f{1}).(subf{1}) = [];
          end
       end
    end    
    
    % merge current subject data into larger struct
    for f = fields' 
        subjDataAll.(f{1}) = mergestruct(subjDataAll.(f{1}),subjData.res.(f{1}));
    end
    
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



function sout = mergestruct(struct1,struct2)

if isempty(struct1) & ~isempty(struct2)
    sout = struct2;
    return
elseif ~isempty(struct1) & isempty(struct2)
    sout = struct1;
    return
end

fields1 = fieldnames(struct1);
fields2 = fieldnames(struct2);

if isequal(fields1,fields2)
    sout = cell2struct(fields1,fields1,1);
    for f = 1:length(fields1)
        sout.(fields1{f}) = horzcat(struct1.(fields1{f}),struct2.(fields1{f}));       
    end
end


