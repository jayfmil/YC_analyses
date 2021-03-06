function aucByCorrectThresh(subjs,params,saveDir,numIters,percentiles)

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_TH1');
end

% load parameters, if not given
if ~exist('params','var') || isempty(params)
    params = TH_multiParams;
end

% number of shuffles
if ~exist('numIters','var') || isempty(numIters)
    numIters = 100;
end

% number of shuffles
if ~exist('percentiles','var') || isempty(percentiles)
    percentiles = 10:2:90;
end

% make base save directory
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

subjData = cell(1,length(subjs));
% for each subject, run classifier for a given correct/incrorrect threshold
for s = 1:length(subjs)
    subj = subjs{s};
    subjFile = fullfile(saveDir,[subj '_aucByThresh.mat']);
    if exist(subjFile,'file')
        subjData{s} = load(subjFile);
    else
        % load events to figure out the threhsolds
        events = get_sub_events('RAM_TH1',subjs{s});
        distErrs = [events(strcmp({events.type},'CHEST')).distErr];
        distErrs_sort = sort(distErrs(~isnan(distErrs)));
%         thresholds = round(distErrs_sort(5)+1):round(distErrs_sort(end-5)-1);
        thresholds = prctile(distErrs,percentiles);
        
        % aucs by threshold, aucs for each permutation by threshold, pvalue at
        % each threshold
        aucs     = NaN(1,length(thresholds));
        aucsPerm = NaN(numIters,length(thresholds));
        p        = NaN(1,length(thresholds));
        for thresh = 1:length(thresholds)
            
            % set the correct threshold to classify
            params.correctThresh = thresholds(thresh);
            params.loadPower  = 0;
            params.doPermute = 0;
            params.saveOutput = 1;
            params.cvField = 'session';
            params.timeBins = [55 120];
            params.freqBins = [];
            params.savePower = 1;
            params.powerPath = '/scratch/jfm2/power16freqs';
            params.basePath = saveDir;
            params.overwrite = 1;
            Cs = logspace(log10(1e-6),log10(1e4),22);
            params.Cs = Cs(7);
            
            subSaveDir = fullfile(saveDir,['thresh_',num2str(params.correctThresh)]);
            subSaveDir = fullfile(saveDir,subj);
            if thresh == 1
                params.savePower = 1;
                params.loadPower = 0;
            else
                params.savePower = 0;
                params.loadPower = 1;
            end
            if ~exist(subSaveDir,'dir');mkdir(subSaveDir);end
            aucs(thresh) = TH1_refactor_phase(subjs{s},params,subSaveDir);
            
            params.loadPower = 1;
            params.doPermute = 1;
            params.saveOutput = 0;
            params.savePower = 0;
            parfor iter = 1:numIters
                aucsPerm(iter,thresh) = TH1_refactor_phase(subjs{s},params,subSaveDir);
            end
            p(thresh) = mean(aucs(thresh) < aucsPerm(:,thresh));
        end
        save(subjFile,'aucs','aucsPerm','p','thresholds','params','distErrs','subj','percentiles');
    end    
    subjData{s} = load(subjFile);
end
keyboard

% directory to save figures
figDir = fullfile(saveDir,'reports');
if ~exist(figDir,'dir');mkdir(figDir);end

% loop over each subject. plot auc as a function of distance threshold and
% keep track of best threshold
median_err = NaN(1,length(subjData));
best_dist  = NaN(1,length(subjData));
ps         = NaN(1,length(subjData));
for s = 1:length(subjData)
    figure(1)
    clf
    
    % plot auc
    plot(subjData{s}.thresholds,subjData{s}.aucs,'linewidth',4)
    
    % create x labels with the thresholds
    pCorr_by_thresh = NaN(1,length(subjData{s}.thresholds));
    xlabels = cell(1,length(subjData{s}.thresholds));
    for t = 1:length(pCorr_by_thresh)      
        pCorr_by_thresh(t) = sum(subjData{s}.distErrs < subjData{s}.thresholds(t))/sum(~isnan(subjData{s}.distErrs));
        xlabels{t} = sprintf('%.1f (%.1f)',subjData{s}.thresholds(t),pCorr_by_thresh(t)*100);
    end
    
    % format figure
    set(gca,'fontsize',20)
    set(gca,'xtick',subjData{s}.thresholds(1:4:t))
    set(gca,'xticklabel',xlabels(1:4:t))
    grid on
    set(gca,'gridlinestyle',':')
    ylabel('AUC','fontsize',24);
    xlabel('Distance Threshold (units)','fontsize',24)
    set(gca,'xticklabelRotation',360-45)
    
    % plot horizontal line at .5
    hold on
    set(gca,'xlim',[subjData{s}.thresholds(1) - 1,subjData{s}.thresholds(end)+1])
    xlim = get(gca,'xlim');
    plot(xlim,[.5 .5],'--k','linewidth',2)
    set(gca,'xlim',xlim)
    
    % plot sig. indicators
    pPos = subjData{s}.p < .025;
    y = max(subjData{s}.aucs) + .025;
    if any(pPos)
        plot(subjData{s}.thresholds(find(pPos)),y,'ok','markerfacecolor','r','markersize',8)
    end
    pNeg =  subjData{s}.p > .975;
    y = min(subjData{s}.aucs) - .025;
    if any(pNeg)
        plot(subjData{s}.thresholds(find(pNeg)),y,'ok','markerfacecolor','b','markersize',8)
    end
    
    % plot vertical line at 13, the default threshold
    ylim = get(gca,'ylim');
    plot([13 13],ylim,'--k','linewidth',2)
    set(gca,'ylim',ylim)
    
    % keep track of best distance, as well as median perofrmance
    best_dist(s) = subjData{s}.thresholds(subjData{s}.aucs==max(subjData{s}.aucs));
    median_err(s) = nanmean(subjData{s}.distErrs);
    ps(s) = subjData{s}.p(subjData{s}.aucs==max(subjData{s}.aucs));
    
    % print
    set(gcf,'paperpositionmode','auto')
    fname = fullfile(figDir,[subjData{s}.subj '_aucByThresh']);
    print('-depsc2','-loose',fname);
end

% scatter of best distance threshold vs median performance
figure(2)
clf
scatter(best_dist, median_err,100,'filled');
[r,p] = corr(best_dist',median_err');
stats = regstats(median_err,best_dist);
hold on
x = [min(best_dist) max(best_dist)];
h = plot(x,x*stats.beta(2)+stats.beta(1),'-k','linewidth',3);
titleStr = sprintf('rho=%.2f, p=%.3f',r,p);
xlabel('Best Classifier Distance (VR units)','fontsize',20);
ylabel('Median Behavioral Error (VR units)','fontsize',20);
grid on
set(gca,'gridlinestyle',':')
title(titleStr)
set(gca,'fontsize',20)
fname = fullfile(figDir,'scatter_classDistByBehError');
print('-depsc2','-loose',fname);
keyboard




























