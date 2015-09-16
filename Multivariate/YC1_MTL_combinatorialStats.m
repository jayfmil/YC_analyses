function YC1_MTL_combinatorialStats(basePath)
% function YC1_MTL_combinatorialStats(basePath)
%
% Run this after running YC1_runMulti_MTL_subregion_wrapper(). Fill this
% in.

% tex directory
saveDir = fullfile(basePath,'reports');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% figure directory
figDir = fullfile(saveDir,'figs');
if ~exist('figDir','dir')
    mkdir(figDir)
end

% all combinations of regions
regions = {'CA1','CA3','DG','EC','PHC','PRC','Sub'};
regions = {'mtl','frontal','parietal','occipital','temporal','limbic'};
combs   = nchoosek(1:length(regions),2);

% get subjects
subjs = get_subs('RAM_YC1');

% initialize matrix to hold all MTL AUC (# subjs x # 1)
mtlAUC = NaN(length(subjs),1);

% initialize matrix to hold single region AUC (# subjs x # regions)
single_regionAUC = NaN(length(subjs),length(regions));

% initialize matrix to hold multiregion region AUC (# subjs x # combinations)
multi_regionAUC = NaN(length(subjs),size(combs,1));

% initialize matrix to hold difference between single region and
% multiregion AUC (# regions x # regions x # subjects)
single_multi_diffAUC = NaN(length(regions),length(regions),length(subjs));

% loop over subjects
for s = 1:length(subjs)
    
    % load the MTL AUC
    regionPath = fullfile(basePath,'MTL');
    regionPath = fullfile(basePath,'all');
    fname      = fullfile(regionPath,'OrigPower',[subjs{s} '_lasso.mat']);
    if exist(fname,'file')
        res = load(fname);
        mtlAUC(s) = res.AUC;
    end
    
    % then loop over each region, store AUC
    for r = 1:length(regions)
        regionPath = fullfile(basePath,regions{r});
        fname      = fullfile(regionPath,'OrigPower',[subjs{s} '_lasso.mat']);
        if exist(fname,'file')
            res = load(fname);
            single_regionAUC(s,r) = res.AUC;
        end
    end
    
    % then loop over all combinations of regions, store AUC
    for c = 1:size(combs,1)
        regionPath = fullfile(basePath,sprintf('%s-%s',regions{combs(c,1)},regions{combs(c,2)}));                
        fname      = fullfile(regionPath,'OrigPower',[subjs{s} '_lasso.mat']);
        if exist(fname,'file')
            res = load(fname);
            multi_regionAUC(s,c) = res.AUC;
        end
    end
    
    % now compute the difference bwteen single and multi region AUC
    for r1 = 1:length(regions)    
        
        % initial region AUC
        auc_single   = single_regionAUC(s,r1);
        otherRegions = setdiff(1:length(regions),r1);
        
        % loop over other regions, exlcuding current
        for r2 = otherRegions
            
            % find where the current region was paired with the other
            % region
            idx = sum(combs == r1 | combs == r2,2)==2;
            
            % find paired region AUC and subtract initial
            auc_multi = multi_regionAUC(s,idx);
            auc_diff  = auc_multi - auc_single;
            
            % store in 3d matrix
            single_multi_diffAUC(r1,r2,s) = auc_diff;            
        end
    end
    
end
keyboard
figs = [];
%% FIGURE 1a - heatmap of AUC differences
% myfig = figure('Position',[1,20,1000,1000]);
clf
a = axes('position',[.1 .5 .35 .45]);
fname = fullfile(figDir,'heatmap');
figs.heatmap = fname;

plotData = nanmean(single_multi_diffAUC,3);
plotData(isnan(plotData)) = -.16;
imagesc(plotData);
c=colormap('jet');
colormap([0 0 0;c])
set(gca,'clim',[-.15 .15]);
xlabel('Added ROI','fontsize',16)
ylabel('Initial ROI','fontsize',16)
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'Layer','bottom');
set(gca,'fontsize',16)
axis square
box on
hc = colorbar('NorthOutside');
xlabel(hc,'\Delta AUC')
set(hc,'fontsize',16)
drawnow
a1=gca;

%% FIGURE 1b - ACI BAR
pos = get(a1,'Position');
a=axes('position',[pos(1)+pos(3)+.1 pos(2) pos(3) pos(4)]);
ACI = squeeze(nanmean(single_multi_diffAUC,2))';
[h,p,c,s] = ttest(ACI);
y = nanmean(ACI);
e = nanstd(ACI)./sqrt(sum(~isnan(ACI)));
x = 1:length(y);

if any(h)
    hold on
    ha=bar(x(h==0),y(h==0),'linewidth',2);
    set(ha,'facecolor',[.7 .7 .7]);
    if any(h & s.tstat>0)
        sig = h & s.tstat>0;
        ha=bar(x(sig),y(sig),'linewidth',2);
        set(ha,'facecolor',[.7 .2 .2]);
    end
    if any(h & s.tstat<0)
        sig = h & s.tstat<0;
        ha=bar(x(sig),y(sig),'linewidth',2);
        set(ha,'facecolor',[.2 .2 .7]);        
    end
else
    hold on
    ha=bar(x,y,'linewidth',2);
    set(ha,'facecolor',[.7 .7 .7]);
end
errorbar(x,y,e,'k','linewidth',2,'linestyle','none');
set(gca,'xtick',1:(length(x)+.5))
set(gca,'xticklabel',regions)
set(gca,'xlim',[0.5 length(x)+.5])
set(gca,'ylim',[-.1 .1])
set(gca,'ytick',[-.1:.05:.1])
set(gca,'XAxisLocation','Top')
set(gca,'YAxisLocation','Right')
ylabel('ACI (\Delta AUC)','fontsize',16)
view([90 90])
grid on
set(gca,'gridlinestyle',':')
set(gca,'fontsize',16)
axis square


%% FIGURE 1c - UCP BAR
pos=get(a1,'Position');
a=axes('position',[pos(1) .1 pos(3) pos(4)]);
UCP = squeeze(nanmean(single_multi_diffAUC,1))';
[h,p,c,s] = ttest(UCP);
y = nanmean(UCP);
e = nanstd(UCP)./sqrt(sum(~isnan(UCP)));
x = 1:length(y);

if any(h)
    hold on
    ha=bar(x(h==0),y(h==0),'linewidth',2);
    set(ha,'facecolor',[.7 .7 .7]);
    if any(h & s.tstat>0)
        sig = h & s.tstat>0;
        ha=bar(x(sig),y(sig),'linewidth',2);
        set(ha,'facecolor',[.7 .2 .2]);
    end
    if any(h & s.tstat<0)
        sig = h & s.tstat<0;
        ha=bar(x(sig),y(sig),'linewidth',2);
        set(ha,'facecolor',[.2 .2 .7]);        
    end
else
    hold on
    ha=bar(x,y,'linewidth',2);
    set(ha,'facecolor',[.7 .7 .7]);
end
errorbar(x,y,e,'k','linewidth',2,'linestyle','none');
set(gca,'xtick',1:(length(x)+.5))
set(gca,'xticklabel',regions)
set(gca,'XTickLabelRotation',45)
set(gca,'xlim',[0.5 length(x)+.5])
set(gca,'ylim',[-.1 .1])
set(gca,'ytick',[-.1:.05:.1])
ylabel('UCP (\Delta AUC)','fontsize',16)
grid on
set(gca,'gridlinestyle',':')
set(gca,'fontsize',16)
axis square


%% FIGURE 1d - PERFORMANCE BY NUMBER OF REGIONS BAR
a=axes('position',[pos(1)+pos(3)+.1 .1 pos(3) pos(4)]);
% fname = fullfile(figDir,'perf');
% figs.perf = fname;

% one region average for each subject
subj_single_auc = nanmean(single_regionAUC,2);

% two region average for each subject
subj_multi_auc = nanmean(multi_regionAUC,2);

% one region, two regions, all mtl
aucs = [subj_single_auc subj_multi_auc mtlAUC];

[h,p,c,s] = ttest(aucs,.5);
y = nanmean(aucs);
e = nanstd(aucs)./sqrt(sum(~isnan(aucs)));
x = 1:length(y);

if any(h)
    hold on
    ha=bar(x(h==0),y(h==0),'linewidth',2);
    set(ha,'facecolor',[.7 .7 .7]);
    if any(h & s.tstat>0)
        sig = h & s.tstat>0;
        ha=bar(x(sig),y(sig),'linewidth',2);
        set(ha,'facecolor',[.7 .2 .2]);
    end
    if any(h & s.tstat<0)
        sig = h & s.tstat<0;
        ha=bar(x(sig),y(sig),'linewidth',2);
        set(ha,'facecolor',[.2 .2 .7]);        
    end
else
    hold on
    ha=bar(x,y,'linewidth',2);
    set(ha,'facecolor',[.7 .7 .7]);
end
errorbar(x,y,e,'k','linewidth',2,'linestyle','none');
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'One','Two','All'})
set(gca,'xlim',[0.5 3.5])
set(gca,'ylim',[.4 .6])
set(gca,'ytick',[.4:.05:.6])
ylabel('AUC','fontsize',16)
xlabel('Number of regions','fontsize',16)
grid on
set(gca,'gridlinestyle',':')
set(gca,'fontsize',16)
axis square
box on
set(gcf,'paperpositionmode','auto') 
% print('-dpng','-loose',fname);
% print('-depsc2','-tiff','-loose',fname);
keyboard
%% FIGURE 2 - Subject counts per region
myfig2 = figure('Position',[1,20,1000,400]);
fname = fullfile(figDir,'counts');
figs.counts = fname;
keyboard
subj_counts = sum(~isnan(single_regionAUC));
ha=bar(1:length(subj_counts),subj_counts,'linewidth',2);
set(ha,'facecolor',[.7 .7 .7]);
ylabel('N (patients)','fontsize',16)
set(gca,'xtick',1:7)
set(gca,'xticklabel',regions)
set(gca,'xlim',[0.5 length(subj_counts)+.5])
grid on
set(gca,'gridlinestyle',':')
set(gca,'fontsize',16)
set(gcf,'paperpositionmode','auto') 
% print('-depsc2','-tiff','-loose',fname);
keyboard


texName = 'YC1_combinatorial_MTL_report.tex';
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
fprintf(fid,'\\rhead{YC1 MTL Combinatorial Report All Electrodes -  Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

fprintf(fid,'\\begin{figure}[!h]\n');
fprintf(fid,'\\centering\n');
fprintf(fid,'\\includegraphics[width=1\\textwidth]{%s}\n',figs.heatmap);
fprintf(fid,'\\includegraphics[width=1\\textwidth]{%s}\n',figs.counts);
% fprintf(fid,'\\caption{%d Subjects: Subject average quartile by time bin. YC1 model from each time bin is applied to the \\textbf{same} time bin in YC2.}\n\n',figs.N(1,1));
fprintf(fid,'\\end{figure}\n\n\n');
fprintf(fid,'\\clearpage\n\n\n');





fprintf(fid,'\\end{document}\n\n\n');







