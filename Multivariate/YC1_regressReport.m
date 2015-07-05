function YC1_regressReport(subjs,dataDir,saveDir)

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

if ~exist('dataDir','var') || isempty(dataDir)
    params = multiParams();
    dataDir = params.regressDir;
end

if ~exist('saveDir','var') || isempty(saveDir)
    saveDir = fullfile(dataDir,'reports');    
end

if ~exist('saveDir','dir')
    mkdir(saveDir);
end

% process each subject
for s = 1:length(subjs);
    fprintf('Processing %s.\n',subjs{s});
    process_subj(subjs{s},dataDir,saveDir)
end

function process_subj(subj,dataDir,saveDir)

% does data exist?
if ~exist(fullfile(dataDir,subj),'dir');
    fprintf('No regression corrected power values for %s.\n',subj)
    return
end

% load tal structure
try
    tal = getBipolarSubjElecs(subj,1,1);
catch
    fprintf('Could not load electrode locations for %s. Aborting.\n',subj)
    return
end

% right now we only care if there are any MTL electrodes
if ~isfield(tal,'locTag') || ~any(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
    fprintf('No MTL electrodes for %s.\n',subj)
    return
end

% make sure all the save directories exist
subjReportDir = fullfile(saveDir,subj);
if ~exist('subjReportDir','dir')
    mkdir(subjReportDir);
end

subjFigDir = fullfile(subjReportDir,'figs');
if ~exist('subjFigDir','dir')
    mkdir(subjFigDir);
end

% load config for freqs and time bins
config = RAM_config('RAM_YC1');

% Setting time bins for convenience:
tEnds = (config.distributedParams.timeWin:...
    config.distributedParams.timeStep:...
    config.postMS+config.priorMS)-config.priorMS;
tStarts = tEnds - config.distributedParams.timeWin + 1;
config.distributedParams.timeBins = [tStarts' tEnds'];

% hipp electrodes
hipp_tal = tal(~cellfun('isempty',regexpi({tal.locTag},['CA1|CA3|DG|sub'])));

% EC electrodes
ec_tal = tal(~cellfun('isempty',regexpi({tal.locTag},['EC'])));

% combined
mtl_tal = [hipp_tal ec_tal];

rois = {'hipp','ec','mtl'};
for roi = rois
    
    % filter electrodes by current region
    switch roi{1}
        case 'hipp'
            region_tal = hipp_tal;
            if isempty(region_tal)
                fprintf('No hipp electrodes for %s.\n',subj)
                continue
            end
        case 'ec'
            region_tal = ec_tal;
            if isempty(region_tal)
                fprintf('No ec electrodes for %s.\n',subj)
                continue
            end
            region_tal = mtl_tal;
            if isempty(region_tal)
                fprintf('No mtl electrodes for %s.\n',subj)
                continue
            end
    end
    
    figs = [];
    for elec = 1:length(region_tal)
               
        % load regression corrected data for current electrode
        elecNum  = tal(elec).channel;
        fname    = sprintf('%s_elec_%d-%d_regStatistics.mat',subj,elecNum(1),elecNum(2));
        eFile    = fullfile(dataDir,subj,fname);        
        elecData = load(eFile);
        figName  = sprintf('%s_elec_%d-%d_',subj,elecNum(1),elecNum(2));
        
        % will hold average data for all electrodes
        if elec == 1
            [pval1_all,pval2_all,pval3_all,...
             tstat1_all,tstat2_all,tstat3_all] ...
             = deal(NaN(length(region_tal),size(elecData.pval1,1),size(elecData.pval1,2)));
        end        
        
        % plot predictor 1 (difficulty) pval
        pval1_all(elec,:,:) = elecData.pval1;
        h=plot_time_by_freq(elecData.pval1,1,config,'Difficulty');
        figs(elec).pval1 = fullfile(subjFigDir,[figName 'pval1']);
        print(figs(elec).pval1,'-loose','-depsc2');
        
        % plot predictor 1 (difficulty) tstat
        tstat1_all(elec,:,:) = elecData.tstat1;
        h=plot_time_by_freq(elecData.tstat1,0,config,'Difficulty');
        figs(elec).tstat1 = fullfile(subjFigDir,[figName 'tstat1']);
        print(figs(elec).tstat1,'-loose','-depsc2');
        
        % plot predictor 2 (learning trial 1 vs 2) pval
        pval2_all(elec,:,:) = elecData.pval2;
        h=plot_time_by_freq(elecData.pval2,1,config,'Learning Trial');
        figs(elec).pval2 = fullfile(subjFigDir,[figName 'pval2']);
        print(figs(elec).pval2,'-loose','-depsc2');
        
        % plot predictor 2 (learning trial 1 vs 2) tstat
        tstat2_all(elec,:,:) = elecData.tstat2;
        h=plot_time_by_freq(elecData.tstat2,0,config,'Learning Trial');
        figs(elec).tstat2 = fullfile(subjFigDir,[figName 'tstat2']);
        print(figs(elec).tstat2,'-loose','-depsc2');        
        
        % plot predictor 3 (trial num) pval
        pval3_all(elec,:,:) = elecData.pval3;
        h=plot_time_by_freq(elecData.pval3,1,config,'Trial Number');
        figs(elec).pval3 = fullfile(subjFigDir,[figName 'pval3']);
        print(figs(elec).pval3,'-loose','-depsc2');
        
        % plot predictor 3 (trial num) tstat
        tstat3_all(elec,:,:) = elecData.tstat3;
        h=plot_time_by_freq(elecData.tstat3,0,config,'Trial Number');
        figs(elec).tstat3 = fullfile(subjFigDir,[figName 'tstat3']);
        print(figs(elec).tstat3,'-loose','-depsc2')        
                
    end
    
    % compute average across all electrode in region
    
    keyboard
    
    
end

function h=plot_time_by_freq(data,isPval,config,titleStr)
clf
if isPval
    imagesc(-log10(data)');axis xy;colormap jet;
    gca.CLim = [-max(abs(h.CLim)) max(abs(h.CLim))];
    h=colorbar;
    h.Label.String = '-log10(p)';
    h.Label.FontSize = 14;
else
    imagesc(data');axis xy;colormap jet;
    h=colorbar;
    h.Label.String = 'tstat';
    h.Label.FontSize = 14;
end

% y info
freqs = config.distributedParams.freQ;
h = gca;
h.YTick      = 1:5:length(freqs);
h.YTickLabel = round(freqs(1:5:end));
ylabel('Frequency (Hz)','fontsize',16)

% x info
binsPerS = 1000/config.distributedParams.timeStep;
bins     = 0:binsPerS:size(data,1);
bins(1)  = 1;
times    = [config.distributedParams.timeBins(bins(1:end-1)+1,1); config.distributedParams.timeBins(end)];
h.XTick  = bins;
h.XTickLabel = round(times/1000);
xlabel('Time (s)','fontsize',16)

h.FontSize = 16;
title(titleStr,'fontsize',14);

% 
% 
% figs = [];
% nTimes = length(subjData.res);
% for t = 1:nTimes
%     
%     % for each time bin, get timebin specific info
%     time   = subjData.res(t).timeBin;
%     yHat   = vertcat(subjData.res(t).yHat{:});
%     mse    = vertcat(subjData.res(t).mse{:});
%     lambda = subjData.res(t).lambda;
%     df     = subjData.res(t).df;
%     
%     % mean yHat for each object position (bc two trials per object)
%     yHat_meanByObj = grpstats(abs(abs(subjData.Y)-abs(yHat)),objGroups);
%     
%     % Make Figure
%     close all
%     figure('units','normalized','paperpositionmode','auto','position',[.2  .2  .5  .25]);
%     titleStr = sprintf('%s: %s, Time = %d-%d, DF = %d, Lambda = %.4f',strrep(subj,'_',' '),ana_name,time(1),time(2),df,lambda);
%     title(titleStr,'fontsize',14)
%     axis off
%         
%     % Y x Prediction
%     ax = axes('Position',[.08 .2 .25 .7]);
%     h=scatter(ax,subjData.Y,yHat);
%     ax.YLim = [0 1];
%     ax.YTick = 0:.2:1;
%     ax.XTick = 0:.2:1;
%     grid on
%     ax.GridLineStyle = ':';
%     ax.GridColor = 'k';
%     xlabel('Performance Factor','fontsize',22)
%     ylabel('Prediction','fontsize',22) 
%     ax.Clipping = 'off';
%     ax.FontSize = 22;
%     
%     % Y x Residuals        
%     ax = axes('Position',[.41 .2 .25 .7]);
%     h = scatter(ax,subjData.Y,subjData.Y-yHat);
%     ax.XTick = 0:.2:1;
%     grid on
%     ax.GridLineStyle = ':';
%     ax.GridColor = 'k';
%     xlabel('Performance Factor','fontsize',22)
%     ylabel('Residual','fontsize',22) 
%     ax.Clipping = 'off';
%     ax.FontSize = 22;
%     
%     % 2D scatter prediction error
%     ax = axes('Position',[.68 .56 .3 .4]);
%     scatter(ax,uniqueObjPos(:,1),uniqueObjPos(:,2),100,yHat_meanByObj,'filled');
%     ax.YLim = [-18 18];
%     ax.XLim = [-32 32];   
%     ax.XTick = [];
%     ax.YTick = [];
%     ax.Clipping = 'off';
%     title('Pred. Error')
%     h = colorbar;
%     h.Location = 'northoutside';    
%     box on
%     
%     % 2D scatter mse
%     ax = axes('Position',[.68 .14 .3 .4]);
%     scatter(ax,uniqueObjPos(:,1),uniqueObjPos(:,2),100,mse,'filled');
%     ax.YLim = [-18 18];
%     ax.XLim = [-32 32];    
%     ax.XTick = [];
%     ax.YTick = [];
%     ax.Clipping = 'off';
%     title('MSE')
%     box on
%     h = colorbar;
%     h.Location = 'northoutside';
%     
%     figDir = fullfile(saveDir,strrep(ana_name,' ','_'));
%     if ~exist(figDir,'dir')
%         mkdir(figDir);
%     end
%     figName = sprintf('%s_timebin%d',subj,t);    
%     figs(t).fname = fullfile(figDir,figName);
%     print(figs(t).fname,'-depsc2','-loose')
% end
% 
% 
% texName = [subj '_' strrep(ana_name,' ','_') '_report.tex'];
% write_texfile(figDir,texName,subj,figs)
% 
% 
% curr_dir = pwd;
% cd(figDir);
% fprintf('Compiling pdf...\n');
% unix(['pdflatex -shell-escape ' fullfile(figDir, texName)]);
% unix(['rm ' texName(1:end-3) 'aux']);
% unix(['rm ' texName(1:end-3) 'log']);
% fprintf('Done!\n');
% cd(curr_dir);


% Start making the tex file
function write_texfile(saveDir,texName, subj, figs)

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
fprintf(fid,'\\lhead{Report: %s }\n',strrep(subj,'_','\_'));
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
for f = 1:length(figs)
    
    if f == 1
        fprintf(fid,'\\begin{figure}[!h]\n');
        fprintf(fid,'\\centering\n');
    end
    fprintf(fid,'\\includegraphics[width=0.99\\textwidth]{%s}\n',figs(f).fname);
    fprintf(fid,'\\vspace{.2in}\n')
    if f == length(figs)
%         fprintf(fid,'\\caption{\\textbf{test}.}\n');
        fprintf(fid,'\\end{figure}\n\n\n');
    
    end
end

fprintf(fid,'\\end{document}\n\n\n');










