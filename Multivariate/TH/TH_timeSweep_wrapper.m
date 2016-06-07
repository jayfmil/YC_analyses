function TH_timeSweep_wrapper(subjs)
%
% For each subject, perform classification, sweeping over the encoding
% period. Permform for power, phase, power and phase, phase lag


% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_TH1');
end

% base save directory
baseDir = '/scratch/jfm2/TH1/multi/timeSweep_8Freqs_badChans';
subDir  = fullfile(baseDir,'t_%d_featType_%d');

% time bins correspond to -500 ms to +1500 ms
Ts = 25:125;
nTimes = length(Ts);

% 0 is power, 1 is phase, 2 is power and phase, 3 is phase lag pairs
feats = 0:3;

% number of permutations at each time point
iters = 100;

% hold aucs and pvals for each subj, time, features
aucsAll = NaN(length(subjs),nTimes,4);
psAll   = NaN(length(subjs),nTimes,4);

totalCalcs = length(Ts)*4;
for s = 1:length(subjs)
    subj     = subjs{s};
    subjFile = fullfile(baseDir,[subj '_gridRes.mat']);
    if exist(subjFile,'file')
        subjData = load(subjFile);
        aucsAll(s,:,:) = subjData.aucsSubj;
        psAll(s,:,:)   = subjData.psSubj;
    else
        fprintf('Processing %s.\n',subjs{s})
        aucsSubj = NaN(1,length(Ts),4);
        psSubj   = NaN(1,length(Ts),4);
        parfor ind = 1:totalCalcs
            [t,feat] = ind2sub([nTimes,4],ind);
            
            params = TH_multiParams();
            params.Cs = params.Cs(7);
            params.powerPath = '/scratch/jfm2/power8freqs';
            %powerPaths = {'/scratch/jfm2/power8freqs','/scratch/jfm2/power4bins_hilbert'};
            params.normType = 'L2';
            params.freqBins = [];
            params.saveOutput = 1;
            params.overwrite = 1;
            params.timeBins = [Ts(t) Ts(t)]
            
            % NEW: EXCLUDE BAD CHANNELS
            params.excludeBadChansWatrous = 1;            
            
            if feat==2
                params.timeBins = Ts(t);
            end
            params.usePhase = feats(feat);
            
            % default to session level cv. Will fall back on trial if only one
            % session.
            params.cvField = 'session';
            
            % save directory
            saveDir = sprintf(subDir,t,feat);
            if ~exist(saveDir,'dir')
                mkdir(saveDir)
            end
            params.basePath = saveDir;
            aucsSubj(ind) = TH1_refactor_phase(subj,params,saveDir);
            
            % permute and recalc auc
            aucsPerm = NaN(1,iters);
            params.loadPower = 1;
            params.doPermute = 1;
            params.overwrite = 0;
            params.saveOutput = 0;
            for iter = 1:iters
                aucsPerm(iter) = TH1_refactor_phase(subj,params,saveDir);
            end
            psSubj(ind) = mean(aucsSubj(ind) < aucsPerm)
        end
        aucsAll(s,:,:) = aucsSubj;
        psAll(s,:,:)   = psSubj;
        save(subjFile,'aucsSubj','psSubj','Ts');
    end
end

fname = fullfile(baseDir,'all_gridRes.mat');
save(fname,'aucsAll','psAll','Ts','subjs');
keyboard


% make report
pos = [0.2    0.2    0.9    0.8];
figure('units','normalized','position',pos);
figs = [];
titles = {'Power','Phase','Power and Phase','Pairwise Phase Lag'};
goodSubjs = find(~isnan(aucsAll(:,1,1)));
figDir = fullfile(baseDir,'report');
if ~exist(figDir,'dir')
    mkdir(figDir)
end
count = 1;
for s = goodSubjs'
    
    ylims = [];
    clf
    for subp = 1:4
        subplot(4,1,subp)
        plot(Ts,squeeze(aucsAll(s,:,subp)),'-k','linewidth',3)
        hold on
        plot(Ts(psAll(s,:,subp)<.05),aucsAll(s,psAll(s,:,subp)<.05,subp),'.r','markersize',20);
        h=plot(Ts(psAll(s,:,subp)<.025),aucsAll(s,psAll(s,:,subp)<.025,subp),'.r','markersize',30);
        set(h,'Color',[.7 0 0]);
        x = (Ts - 50) * 20;
        set(gca,'xtick',Ts(1:10:101));
        set(gca,'xticklabel',x(1:10:101));
        
        grid on
        set(gca,'gridlinestyle',':')
        if subp==4
            xlabel('Time (ms)','fontsize',16)
        end
        ylabel('AUC','fontsize',20)
        set(gca,'fontsize',20)
        xlim = get(gca,'xlim');
        plot(xlim,[.5 .5],'--k','linewidth',2)
        set(gca,'xlim',[20 130])
        
        ylim = get(gca,'ylim');
        ylims = vertcat(ylims,ylim);
        plot([Ts(x==0) Ts(x==0)],[0 1],'--k','linewidth',2)
        set(gca,'ylim',ylim)        
        title(titles{subp})        
                        
    end
    ylims = max(ylims(:)); 
    for subp = 1:4
        subplot(4,1,subp)
%         set(gca,'ylim',[.5-(ylims-.5) ylims])
        set(gca,'ylim',[.2 .8]);
        set(gca,'ytick',[.25 .5 .75]);
    end
    set(gcf,'paperpositionmode','auto')
    fname = fullfile(figDir,[subjs{s} '_auc.eps']);
    print('-depsc2','-loose',fname);
    figs(count).aucFig = fname;
    figs(count).subj = subjs{s};
    count = count+1;
end


texName = 'timeSweep_report.tex';
write_texfile(figDir,texName,figs)

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
    fprintf(fid,'\\includegraphics[width=0.99\\textwidth]{%s}\n',figs(s).aucFig);    
        
    fprintf(fid,'\\caption{\\textbf{%s.}}\n',strrep(figs(s).subj,'_',' '));
    fprintf(fid,'\\end{figure}\n\n\n');
    if mod(s,2) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
end

fprintf(fid,'\\end{document}\n\n\n');



