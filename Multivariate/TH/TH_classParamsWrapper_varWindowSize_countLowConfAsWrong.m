function TH_classParamsWrapper_varWindowSize_countLowConfAsWrong(subjs,just_loso)
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
subjs = subjs(~strcmp(subjs,'TJ088'));

% do monpolar (really, average rereference) or bipolar
doBipol = 1

% base save directory
baseDir = '/scratch/jfm2/TH1/multi/gridsearch_varWindowSize_countLowAsWrong'
if ~doBipol
    baseDir = [baseDir '_avg_reref'];
end
subDir  = fullfile(baseDir,'w_%d_b_%d_f_%d_r_%d');

% time window from 100 ms to 
windowSizes  = 5:10:190;
binCenters = mean([[5:windowSizes(1):195-windowSizes(1)]' [5+windowSizes(1):windowSizes(1):195]'],2);

allTimeBins = cell(length(windowSizes),length(binCenters));
for window = 1:size(allTimeBins,1)
    winHalfLength = windowSizes(window)/2;
    for time = 1:size(allTimeBins,2)
        tStart = binCenters(time) - winHalfLength;
        tEnd   = binCenters(time) + winHalfLength;
        if tStart <= 0 || tEnd > 195
            continue
        else
           allTimeBins{window,time} = [tStart tEnd];
        end
    end
end

Cs = logspace(log10(1e-6),log10(1e4),22);
Cs = Cs(7);

powerPaths = {'/scratch/jfm2/power4bins_hilbert','/scratch/jfm2/power8freqs','/scratch/jfm2/power8freqs_new','/scratch/jfm2/power16freqs','/scratch/jfm2/power50freqs'};
if ~doBipol
    powerPaths = strcat(powerPaths,'_avg_reref');
end
regions    = {'all'};%,'mtl'};

nBinCenters = size(allTimeBins,2);
nWindows    = size(allTimeBins,1);
nFreqs   = length(powerPaths);
nRegions = length(regions);

aucsAll = NaN(nWindows,nBinCenters,nFreqs,nRegions,length(subjs));
totalCalcs = nBinCenters*nWindows*nFreqs*nRegions;
for s = 1:length(subjs)
    subj     = subjs{s};
    subjFile = fullfile(baseDir,[subj '_gridRes.mat']);
    if just_loso
        events = get_sub_events('RAM_TH1',subj);
        if length(unique([events.session])) == 1
            continue
        end
    end    
    if exist(subjFile,'file')
        subjData = load(subjFile);
        aucsAll(:,:,:,:,s) = subjData.aucsSubj;
    else
        fprintf('Processing %s.\n',subjs{s})
        aucsSubj = NaN(nWindows,nBinCenters,nFreqs,nRegions);
        parfor ind = 1:totalCalcs
            
            [w,b,f,r] = ind2sub([nWindows,nBinCenters,nFreqs,nRegions],ind);
            

            timeWin = allTimeBins{w,b};
            if ~isempty(timeWin)
                
                % set parameters for this combination of parameters
                params = TH_multiParams;
                params.timeBins = timeWin;
                params.Cs = Cs;
                params.freqBins = [];
                params.powerPath = powerPaths{f};
                params.region = regions{r};
                params.savePower = 0;
                params.overwrite = 1;
                                
                if ~doBipol
                    params.doBipol = doBipol;              
                end
                
                params.countLowConfAsWrong = 1;
                
                % NEW: EXCLUDE BAD CHANNELS
                %             params.excludeBadChansWatrous = 1;
                
                % default to session level cv. Will fall back on trial if only one
                % session.
                params.cvField = 'session';
                
                % save directory
                saveDir = sprintf(subDir,w,b,f,r);
                if ~exist(saveDir,'dir')
                    mkdir(saveDir)
                end
                
                params.basePath = saveDir;
                aucsSubj(ind) = TH1_refactor_phase(subj,params,saveDir);                
            end
        end
        aucsAll(:,:,:,:,s) = aucsSubj;
        save(subjFile,'aucsSubj');        
    end
end
fname = fullfile(baseDir,'gridRes.mat');
save(fname,'aucsAll','allTimeBins','Cs','windowSizes','binCenters','powerPaths','regions','subjs')

% make some visualizations. 
binCentersMS  = (binCenters - 100) * 20;
windowSizesMS = windowSizes*20;

figs = [];
count = 1;
freqLabels = {'4 bands','8 Frequencies (3-180)','8 Frequencies (1-200)','16 Frequencies (1-200)','50 Frequencies (1-200)'};
figDir = fullfile(baseDir,'report');
if ~exist(figDir,'dir')
    mkdir(figDir)
end

% WTF. Where are the zeros coming from.
aucsAll(aucsAll==0) = NaN;
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
        sanePColor(nanmean(aucs,3));
        set(gca,'xtick',[1:2:nBinCenters])
        set(gca,'xticklabel',binCentersMS(1:2:end))        
        xlabel('Window Center (ms)','fontsize',16)
        set(gca,'xticklabelrotation',45)
        set(gca,'ytick',[1:2:nWindows])
        set(gca,'yticklabel',windowSizesMS(1:2:end))
        ylabel('Window Size (ms)','fontsize',16)
        set(gca,'fontsize',16)
        grid on
        set(gca,'gridlinestyle',':')
        z = find(binCentersMS > 0,1,'first') - .5;
        hold on
        plot([z z],ylim,'-k','linewidth',2);
        h=colorbar;
        aucsMean = nanmean(aucs,3);
        [row,col] = find(aucsMean == max(aucsMean(:)));
        tStart = binCentersMS(col) - windowSizesMS(row)/2;
        tEnd   = binCentersMS(col) + windowSizesMS(row)/2;
        titleStr = sprintf('Max AUC: %.3f, Center = %d ms, size = %d ms\n range = %d - %d ms',max(aucsMean(:)),binCentersMS(col),windowSizesMS(row),tStart,tEnd);
        title(titleStr)        
        fname = fullfile(figDir,sprintf('auc_f%d_r%d.eps',f,r));
        print('-depsc2','-loose',fname);
        figs(count).aucFig = fname;

        % significance map
        figure(2)
        clf
        
        ts_mask = ts;
        ts_mask(ps >= .05) = 0;
        sanePColor(ts_mask);
        clim = get(gca,'clim');
        clim = [-max(abs(clim)) max(abs(clim))];
        if clim(1) < -6
            clim = [-6 6];
        end
        set(gca,'clim',clim);
        set(gca,'xtick',[1:2:nBinCenters])
        set(gca,'xticklabel',binCentersMS(1:2:end))        
        xlabel('Window Center (ms)','fontsize',16)
        set(gca,'xticklabelrotation',45)
        set(gca,'ytick',[1:2:nWindows])
        set(gca,'yticklabel',windowSizesMS(1:2:end))
        ylabel('Window Size (ms)','fontsize',16)
        set(gca,'fontsize',16)
        grid on
        set(gca,'gridlinestyle',':')
        z = find(binCentersMS > 0,1,'first') - .5;
        hold on
        plot([z z],ylim,'-w','linewidth',2);
        h=colorbar;
        title('t-stat')
        colormap('jet');
        cmap = colormap;
        cmap(32:33,:) = [0 0 0;0 0 0];
        colormap(cmap)
        fname = fullfile(figDir,sprintf('tstat_f%d_r%d.eps',f,r));
        print('-depsc2','-loose',fname);
        figs(count).tstatFig = fname;                                  
        
        
        for s = 1:length(subjs)
            figure(1)
            clf
            aucs  = squeeze(aucsAll(:,:,f,r,s));
            sanePColor(aucs);
            set(gca,'xtick',[1:2:nBinCenters])
            set(gca,'xticklabel',binCentersMS(1:2:end))
            xlabel('Window Center (ms)','fontsize',16)
            set(gca,'xticklabelrotation',45)
            set(gca,'ytick',[1:2:nWindows])
            set(gca,'yticklabel',windowSizesMS(1:2:end))
            ylabel('Window Size (ms)','fontsize',16)
            set(gca,'fontsize',16)
            grid on
            set(gca,'gridlinestyle',':')
            z = find(binCentersMS > 0,1,'first') - .5;
            hold on
            plot([z z],ylim,'-k','linewidth',2);
            h=colorbar;            
            [row,col] = find(aucs == max(aucs(:)));
            if ~isempty(row)
                row = row(1);
                col = col(1);
                tStart = binCentersMS(col) - windowSizesMS(row)/2;
                tEnd   = binCentersMS(col) + windowSizesMS(row)/2;
                subjStr = strrep(subjs{s},'_',' ');
                titleStr = sprintf('%s Max AUC: %.3f, Center = %d ms, size = %d ms\n range = %d - %d ms',subjStr,max(aucs(:)),binCentersMS(col),windowSizesMS(row),tStart,tEnd);
                title(titleStr)
            end
            fname = fullfile(figDir,sprintf('auc_f%d_r%d_%s.eps',f,r,subjs{s}));
            print('-depsc2','-loose',fname);
            figs(count).(subjs{s}) = fname;
        end
        
        count = count+1;
    end
end


texName = 'auc_triangle_countLowAsWrong_bipol';
if doBipol == 0
    texName = 'auc_triangle_countLowAsWrong_avg_ref';
end
if just_loso
    texName = [texName '_loso'];
end
texName = [texName '.tex'];
write_texfile(figDir,texName,figs,subjs)

curr_dir = pwd;
cd(figDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(figDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);
keyboard


% Start making the tex file
function write_texfile(saveDir,texName, figs, subjs)

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
    fprintf(fid,'\\includegraphics[width=0.47\\textwidth]{%s}\n',figs(s).aucFig);
    fprintf(fid,'\\includegraphics[width=0.47\\textwidth]{%s}\n',figs(s).tstatFig);
        
    fprintf(fid,'\\caption{\\textbf{%s.} %d patients: %s}\n',figs(s).label,figs(s).nSubj,figs(s).subjs);
    fprintf(fid,'\\end{figure}\n\n\n');
    if mod(s,3) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
end

for f = 1:length(figs)
    
    for s = 1:length(subjs)
        
        if mod(s,2) == 1
            fprintf(fid,'\\begin{figure}[!h]\n');
            fprintf(fid,'\\centering\n');
        end
        
        fprintf(fid,'\\includegraphics[width=0.47\\textwidth]{%s}\n',figs(f).(subjs{s}));
        
        if mod(s,2) == 0 || s == length(subjs)
            fprintf(fid,'\\caption{\\textbf{%s.}}\n',figs(f).label);
            fprintf(fid,'\\end{figure}\n\n\n');            
%             fprintf(fid,'\\clearpage\n\n\n');
        end
    end
    fprintf(fid,'\\clearpage\n\n\n');
end


fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);



















