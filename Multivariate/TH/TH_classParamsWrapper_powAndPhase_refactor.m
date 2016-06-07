function TH_classParamsWrapper_powAndPhase_refactor(subjs)
%
% For each subject, perform classification, iterating over a number of
% paramters:
%
%           time period
%           penalty parameter (C)
%           power (8 freqs vs 4 bands)
%           regions (mtl vs all)
%

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_TH1');
end


% base save directory
baseDir = '/scratch/jfm2/TH1/multi/gridsearch_powAndPhase';
subDir  = fullfile(baseDir,'t_%d_c_%d_f_%d_r_%d');

% possible parameters
timeBins = 25:5:125;
Cs         = logspace(log10(1e-6),log10(1e4),22);
powerPaths = {'/scratch/jfm2/power8freqs','/scratch/jfm2/power4bins_hilbert'};
regions    = {'all','mtl'};

nTimes   = length(timeBins);
nCs      = length(Cs);
nFreqs   = length(powerPaths);
nRegions = length(regions);

aucsAll = NaN(nTimes,nCs,nFreqs,nRegions,length(subjs));
totalCalcs = nTimes*nCs*nFreqs*nRegions;
for s = 1:length(subjs)
    subj     = subjs{s};
    subjFile = fullfile(baseDir,[subj '_gridRes.mat']);
    if exist(subjFile,'file')
        subjData = load(subjFile);
        aucsAll(:,:,:,:,s) = subjData.aucsSubj;
    else
        fprintf('Processing %s.\n',subjs{s})
        aucsSubj = NaN(nTimes,nCs,nFreqs,nRegions);
        parfor ind = 1:totalCalcs
            [t,c,f,r] = ind2sub([nTimes,nCs,nFreqs,nRegions],ind);
            
            % set parameters for this combination of parameters
            params = TH_multiParams;
            params.timeBins = [timeBins(t) timeBins(t)];
            params.Cs = Cs(c);
            params.freqBins = [];
            params.powerPath = powerPaths{f};
            params.region = regions{r};
            params.savePower = 0;
            params.overwrite = 1;
            params.usePhase = 2;
            
            % default to session level cv. Will fall back on trial if only one
            % session.
            params.cvField = 'session';
            
            % save directory
            saveDir = sprintf(subDir,t,c,f,r);
            if ~exist(saveDir,'dir')
                mkdir(saveDir)
            end
            params.basePath = saveDir;
            aucsSubj(ind) = TH1_refactor_phase(subj,params,saveDir);
        end
        aucsAll(:,:,:,:,s) = aucsSubj;
        save(subjFile,'aucsSubj');
    end
end
fname = fullfile(baseDir,'gridRes.mat');
save(fname,'aucsAll','timeBins','Cs','powerPaths','regions','subjs')

% make some visualizations. imagesc time x penalty for each power and each
% region
% fun = @(x) [num2str((x(1)-50)*20),'-',num2str((x(2)-50)*20)];
timeBins_str = cellfun(@num2str,num2cell((timeBins-50).*20),'UniformOutput',false);
figs = [];
count = 1;
freqLabels = {'8 Frequencies','4 bands'};
figDir = fullfile(baseDir,'report');
if ~exist(figDir,'dir')
    mkdir(figDir)
end
for f = 1:nFreqs
    for r = 1:nRegions
        
        % aucs for this space
        aucs  = squeeze(aucsAll(:,:,f,r,:));
        fun2 = @(x) strrep(x,'_',' ');
        subjs_used = subjs(~isnan(squeeze(aucs(1,1,:))));
        subjs_used = cellfun(fun2,subjs_used,'Uniformoutput',false);
        nSubj = length(subjs_used);
        figs(count).label = [freqLabels{f}, ' ', regions{r}, ' electrodes'];
        figs(count).nSubj = nSubj;
        figs(count).subjs = strjoin(subjs_used,', ');
        
        % ttest vs .5
        [h,p,c,s] = ttest(permute(aucs,[3 1 2]),.5);
        ts = squeeze(s.tstat);
        ps = squeeze(p);
        
        % imagesc
        figure(1)
        clf
        imagesc(nanmean(aucs,3))
        set(gca,'ytick',[1:length(timeBins_str)])
        set(gca,'yticklabel',timeBins_str)
        set(gca,'xticklabel',' ')
        xlabel('Regularization Strength (High to Low)','fontsize',16)
        ylabel('Encoding Time (ms)','fontsize',16)
        set(gca,'fontsize',16)
        h=colorbar;
        title('AUC')
        fname = fullfile(figDir,sprintf('auc_f%d_r%d.png',f,r));
        print('-dpng',fname);
        figs(count).aucFig = fname;

        % significance map
        figure(2)
        clf
        
        ts_mask = ts;
        ts_mask(ps >= .05) = 0;
        imagesc(ts_mask)
        clim = get(gca,'clim');
        clim = [-max(abs(clim)) max(abs(clim))];
        set(gca,'clim',clim)
        set(gca,'ytick',[1:length(timeBins_str)])
        set(gca,'yticklabel',timeBins_str)
        set(gca,'xticklabel',' ')
        xlabel('Regularization Strength (High to Low)','fontsize',16)
        ylabel('Encoding Time (ms)','fontsize',16)
        set(gca,'fontsize',16)
        h=colorbar;
        title('t-stat')
        colormap('jet');
        cmap = colormap;
        cmap(32:33,:) = [1 1 1;1 1 1];
        colormap(cmap)
        fname = fullfile(figDir,sprintf('tstat_f%d_r%d.png',f,r));
        print('-dpng',fname);
        figs(count).tstatFig = fname;        
        count = count+1;
    end
end

% do power and phase lag



texName = 'gridReport_phase.tex';
write_texfile(figDir,texName,figs)

curr_dir = pwd;
cd(figDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(figDir, texName)]);
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
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).aucFig);
    fprintf(fid,'\\includegraphics[width=0.4\\textwidth]{%s}\n',figs(s).tstatFig);
        
    fprintf(fid,'\\caption{\\textbf{%s.} %d patients: %s}\n',figs(s).label,figs(s).nSubj,figs(s).subjs);
    fprintf(fid,'\\end{figure}\n\n\n');
    if mod(s,2) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
end

fprintf(fid,'\\end{document}\n\n\n');




















