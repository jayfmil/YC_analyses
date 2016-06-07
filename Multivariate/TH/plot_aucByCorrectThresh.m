function plot_aucByCorrectThresh(subjs,dataDir)


% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_TH1');
end

for s = 1:1;%length(subjs)
    
    subjFile = fullfile(dataDir,[subjs{s} '_aucByThresh']);
    subjData = load(subjFile);
    
    
    
end





clf
ax = axes;
% set(ax,'Position',[0.0700 0.1100 0.5550 0.8150])
plot(subjData.thresholds,subjData.aucs,'linewidth',4)
% set(gca,'ylim',[.25 .75])

pCorr_by_thresh = NaN(1,length(subjData.thresholds));
xlabels = cell(1,length(subjData.thresholds));


for t = 1:length(pCorr_by_thresh)
    
    threshDir = fullfile(dataDir,['thresh_',num2str(subjData.thresholds(t))]);
    threshFile = fullfile(threshDir,[subjs{s} '_class_pow.mat']);
    threshData = load(threshFile);
    
    if t == 1
        pred_by_thresh  = NaN(length(threshData.distanceErrs),length(subjData.thresholds));
    end
    pred_by_thresh(:,t) = vertcat(threshData.res.yPreds{:})==1;
    
    pCorr_by_thresh(t) = mean(vertcat(threshData.res.yTests{:})==1);
    xlabels{t} = sprintf('%d (%.0f)',subjData.thresholds(t),pCorr_by_thresh(t)*100);
    %xlabels{t} = sprintf('%d',subjData.thresholds(t));
end

% add percent high confidence

set(gca,'fontsize',20)
set(gca,'xtick',subjData.thresholds(1:2:t))
set(gca,'xticklabel',xlabels(1:2:t))
grid on
set(gca,'gridlinestyle',':')
ylabel('AUC','fontsize',24);
xlabel('Distance Threshold (units)','fontsize',24)
set(gca,'xticklabelRotation',360-45)

hold on
set(gca,'xlim',[subjData.thresholds(1) - 1,subjData.thresholds(end)+1])
xlim = get(gca,'xlim');
plot(xlim,[.5 .5],'--k','linewidth',2)
set(gca,'xlim',xlim)

pPos = subjData.p < .025;
y = max(subjData.aucs) + .025;
if any(pPos)
    plot(subjData.thresholds(find(pPos)),y,'ok','markerfacecolor','r','markersize',8)
end

pNeg =  subjData.p > .975;
y = min(subjData.aucs) - .025;
if any(pNeg)
    plot(subjData.thresholds(find(pNeg)),y,'ok','markerfacecolor','b','markersize',8)
end
ylim = get(gca,'ylim');
plot([13 13],ylim,'--k','linewidth',2)
set(gca,'ylim',ylim)

% title(subjs{s})
set(gcf,'paperpositionmode','auto')
keyboard

figure(2)
clf
events = get_sub_events('RAM_TH1',subjs{s});
events = events(threshData.eventsToUse);
xs     = [events.locationX];
ys     = [events.locationY];
xs_rec = [events.chosenLocationX];
ys_rec = [events.chosenLocationY];
distErrs = [events.distErr]';

rec = [events.recalled]==1;
scatter(xs(rec),ys(rec))
hold on
scatter(xs(~rec),ys(~rec))

x_offset = xs_rec - xs;
y_offset = ys_rec - ys;

for t = 1:length(subjData.thresholds)
   clf
   plot(x_offset(pred_by_thresh(:,t)==1),y_offset(pred_by_thresh(:,t)==1),'.r','markersize',24)
   hold on
   plot(x_offset(pred_by_thresh(:,t)==0),y_offset(pred_by_thresh(:,t)==0),'.b','markersize',24)
   viscircles([0,0],subjData.thresholds(t));
   set(gca,'xlim',[-70 70])
   set(gca,'ylim',[-70 70]);
   axis square
   
   pCorr_in = mean(pred_by_thresh(distErrs <= subjData.thresholds(t))==1);
   pCorr_out = mean(pred_by_thresh(distErrs > subjData.thresholds(t))==0);
   pCorr = mean(pred_by_thresh(:,t) == (distErrs <= subjData.thresholds(t))==1);
   titleStr = sprintf('Thresh %d: in = %.3f, out = %.3f, all = %.3f',subjData.thresholds(t),pCorr_in,pCorr_out,pCorr);
   title(titleStr)
   pause
end


% ax = axes;
% set(ax,'Position',[0.6800 0.1100 0.280 0.8150])
% plot(thresholds,pCorr_by_thresh,'linewidth',2)
% set(gca,'fontsize',16)
% xlabel('Distance Threshold','fontsize',20);
% ylabel('Percent Correct','fontsize',20);
% grid on
% set(gca,'gridlinestyle',':')
% hold on
% ylim = get(gca,'ylim');
% plot([13 13],ylim,'--k','linewidth',2)
% set(gca,'ylim',ylim)

%print('-depsc2','-loose',subjs{s})




