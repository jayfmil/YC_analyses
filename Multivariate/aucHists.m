function [perfs,aucs,subjs,nTrials,pfRange,pen,toUse] = aucHists(subjs,params)

% if not given, use default params
if ~exist('params','var') || isempty(params)
    params = multiParams();
end

% overwrite exisiting figures?
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

% tex directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
dataDir = fullfile(params.basePath,f(params.useCorrectedPower,y));
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

f = @(x) max(x(:));

perfs  = NaN(length(subjs),1);
aucs  = NaN(length(subjs),1);
pvals = NaN(length(subjs),1);
pfRange = NaN(length(subjs),1);
pfMed = NaN(length(subjs),1);
pfMedDiff = NaN(length(subjs),1);
nTrials = NaN(length(subjs),1);
pen = NaN(length(subjs),1);
ts = NaN(length(subjs),1);
penMean = NaN(length(subjs),1);
trainAucMean = NaN(length(subjs),1);
tMost = NaN(length(subjs),1);

testLocSim1 = NaN(length(subjs),1);
testLocSim2 = NaN(length(subjs),1);

for s = 1:length(subjs)
    subj = subjs{s};
    
    % see if files exist for subject. if not, continue
    if ~isfield(params,'usePhase') || params.usePhase == 1
        lassoFile  = fullfile(dataDir,[subj '_lasso_phase.mat']);
    elseif params.usePhase == 0
        lassoFile  = fullfile(dataDir,[subj '_lasso_pow.mat']);
    elseif params.usePhase == 2
        lassoFile  = fullfile(dataDir,[subj '_lasso_powphase.mat']);
    end
    if ~exist(lassoFile,'file')
        fprintf('Lasso not found for %s.\n',subj)
        continue
    end
    subjData = load(lassoFile);
    if length(subjData.Y) < 48
        continue
    end
    
    events = get_sub_events('RAM_YC1',subj);
    events = addExtraYCFields(events);
    pf = [events(strcmp({events.type},'NAV_TEST')).testError];
    pfRange(s) = std(1-pf);
    pfMed(s) = median(1-pf);
    
    pfMedDiff(s) = abs(median(pf(1-pf > median(1-pf))) - median(pf(1-pf < median(1-pf))));
    
    
    
    perfs(s) = mean(vertcat(subjData.res.err{:}));
    aucs(s)  = subjData.AUC;
    pvals(s) = 1-binocdf(sum(vertcat(subjData.res.err{:})),length(vertcat(subjData.res.err{:})),.5);
    nTrials(s) = length(subjData.Y);
    pen(s) = std(log10([subjData.res.C{:}]));
    penMean(s) = mean(log10([subjData.res.C{:}]));
    
    ts(s) = std(([subjData.res.tBest{:}]));
    
    t=[subjData.res.tBest{:}];
    tMost(s) = max(grpstats(t,t,'numel')/length(t));
    
    trainingAucs = cellfun(f,subjData.res.aucs);
    trainAucMean(s) = mean(trainingAucs);
    
    
    testEvents = find(strcmp({events.type},'NAV_TEST'));
    learnOneTestDist = NaN(1,length(testEvents));
    learnTwoTestDist = NaN(1,length(testEvents));
    for i = 1:length(learnTwoTestDist)
        learnLoc1 = events(testEvents(i)-2).startLocs;
        learnLoc2 = events(testEvents(i)-1).startLocs;
        testLoc  = events(testEvents(i)).startLocs;
        learnOneTestDist(i) = calc_YC_error(learnLoc1,testLoc);
        learnTwoTestDist(i) = calc_YC_error(learnLoc2,testLoc);
    end
    [r,p]=corr(learnOneTestDist',1-pf');
    testLocSim1(s) = r;
    
    [r,p]=corr(learnTwoTestDist',1-pf');
    res = regstats(1-pf,learnTwoTestDist,'linear',{'rsquare'});
    testLocSim2(s) = r;
end
keyboard
% plot initial hist
figure
clf
negSig = perfs(pvals > .975);
posSig = perfs(pvals < .025);
notSig = perfs(pvals >= .025 & pvals <= .975);
n1 = histc(negSig,0.025:.025:.975); if iscolumn(n1);n1=n1';end
n2 = histc(posSig,0.025:.025:.975); if iscolumn(n2);n2=n2';end
n3 = histc(notSig,0.025:.025:.975); if iscolumn(n3);n3=n3';end
h  = bar([.05:.025:1],[n1' n2' n3'],1,'stacked','linewidth',2);
set(h(3),'FaceColor','w');
set(h(1),'FaceColor',[50,124,203]/255);
set(h(2),'FaceColor',[226 55 67]/255);
xlabel('Prob. Correct','fontsize',20)
ylabel('Count','fontsize',20)
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)

% plot mean
m   = nanmean(perfs);
med = nanmedian(perfs);
hold on
ylim = [0 round(max(get(gca,'ylim'))+1)];
plot(.5*[1 1],ylim','-k','linewidth',2)
plot(m*[1 1],ylim','--r','linewidth',2)
set(gca,'ylim',ylim)
set(gca,'xtick',0:.25:1)
set(gca,'xlim',[0 1])
for i = 1:3
    uistack(h(i),'top')
end

% group ttest against .5
[h,p] = ttest(perfs,.5);
titleStr = sprintf('mean: %.2f, median %.2f, p: %.3f',m,med,p);
title(titleStr,'fontweight','normal','fontsize',14)
print('-depsc2','-loose',fullfile(figDir,'perfHist'))

figure
clf
scatter(pen,perfs*100,30)
ylabel('% Correct','fontsize',20)
xlabel('std(log10(Cs))','fontsize',20);
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)
hold on
% s = regstats(pen,perfs);
print('-depsc2','-loose',fullfile(figDir,'pCorrVsPen'))

figure
clf
scatter(penMean,perfs*100,30)
ylabel('% Correct','fontsize',20)
xlabel('mean(log10(Cs))','fontsize',20);
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)
hold on
% s = regstats(pen,perfs);
print('-depsc2','-loose',fullfile(figDir,'pCorrVsMeanPen'))


figure
clf
scatter(ts,perfs*100,30)
ylabel('% Correct','fontsize',20)
xlabel('std(t)','fontsize',20);
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)
hold on
% s = regstats(pen,perfs);
print('-depsc2','-loose',fullfile(figDir,'pCorrVsT'))

figure
clf
scatter(nTrials,perfs*100,30)
ylabel('% Correct','fontsize',20)
xlabel('# Trials','fontsize',20);
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)
hold on
% s = regstats(pen,perfs);
print('-depsc2','-loose',fullfile(figDir,'pCorrVsnTrials'))


figure
clf
scatter(pfRange,perfs*100,30)
ylabel('% Correct','fontsize',20)
xlabel('std(Performance Factor)','fontsize',20);
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)
hold on
% s = regstats(pen,perfs);
print('-depsc2','-loose',fullfile(figDir,'pCorrVsStdPF'))

figure
clf
scatter(pfMed,perfs*100,30)
xlabel('median(Performance Factor)','fontsize',20)
ylabel('% Correct','fontsize',20);
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)
hold on
% s = regstats(pen,perfs);
print('-depsc2','-loose',fullfile(figDir,'pCorrVsMedPF'))


figure
clf
scatter(trainAucMean,perfs*100,30)
xlabel('mean(Training AUC)','fontsize',20)
ylabel('% Correct','fontsize',20);
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)
hold on
% s = regstats(pen,perfs);
print('-depsc2','-loose',fullfile(figDir,'pCorrVsTrainAUC'))


% figure
% clf
% scatter(testLocSim,perfs*100,30)
% xlabel('Learn 2 Test Sim','fontsize',20)
% ylabel('% Correct','fontsize',20);
% grid on
% set(gca,'gridlinestyle',':');
% set(gca,'fontsize',20)
% hold on
% % s = regstats(pen,perfs);
% print('-depsc2','-loose',fullfile(figDir,'pCorrVsTrainAUC'))




figure
clf
scatter(tMost,perfs*100,30)
xlabel('% Mode Timepoint','fontsize',20)
ylabel('% Correct','fontsize',20);
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)
hold on
% s = regstats(pen,perfs);
print('-depsc2','-loose',fullfile(figDir,'pCorrVsTMost'))

% plot hist only high t most
toUse = tMost > .5;
figure
clf

negSig = perfs(pvals > .975 & toUse);
posSig = perfs(pvals < .025 & toUse);
notSig = perfs(pvals >= .025 & pvals <= .975 & toUse);
n1 = histc(negSig,0.025:.025:.975); if iscolumn(n1);n1=n1';end
n2 = histc(posSig,0.025:.025:.975); if iscolumn(n2);n2=n2';end
n3 = histc(notSig,0.025:.025:.975); if iscolumn(n3);n3=n3';end
h  = bar([.05:.025:1],[n1' n2' n3'],1,'stacked','linewidth',2);
set(h(3),'FaceColor','w');
set(h(1),'FaceColor',[50,124,203]/255);
set(h(2),'FaceColor',[226 55 67]/255);
xlabel('Prob. Correct','fontsize',20)
ylabel('Count','fontsize',20)
grid on
set(gca,'gridlinestyle',':');
set(gca,'fontsize',20)

% plot mean
m   = nanmean(perfs(toUse));
med = nanmedian(perfs(toUse));
hold on
ylim = [0 round(max(get(gca,'ylim'))+1)];
plot(.5*[1 1],ylim','-k','linewidth',2)
plot(m*[1 1],ylim','--r','linewidth',2)
set(gca,'ylim',ylim)
set(gca,'xtick',0:.25:1)
set(gca,'xlim',[0 1])
for i = 1:3
    uistack(h(i),'top')
end

% group ttest against .5
[h,p] = ttest(perfs(toUse),.5);
titleStr = sprintf('mean: %.2f, median %.2f, p: %.3f',m,med,p);
title(titleStr,'fontweight','normal','fontsize',14)
print('-depsc2','-loose',fullfile(figDir,'perfHist_onlyHighT'))
keyboard



















