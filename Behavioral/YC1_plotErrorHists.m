function YC1_plotErrorHists()
% function YC1_plotErrorHists()
%
% Makes histograms of distance error (both normalized and euclidean) for
% the entire YC1 dataset. Method 1 (figures 1 and 2): All errors from all
% subjects, no averaging. Method 2 (figures 3 and 4): averaged across
% subjects.

% load all errors
[allErr,allObjs,allEuc,subj] = YC1_loadAllSubjErrors(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% METHOD 1: no averaging %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized distance error
figure(1)
[n,x] = hist(allErr,100);   
n = n/sum(n);
bar(x,n,1,'facecolor',[.5 .5 .5])
grid on
xlabel('Norm. Distance Error','fontsize',20)
ylabel('Probability','fontsize',20)               
set(gca,'fontsize',20)

% euclidean distance error
figure(2)
[n,x] = hist(allEuc,100);   
n = n/sum(n);
bar(x,n,1,'facecolor',[.5 .5 .5])
grid on
xlabel('Euc. Distance Error','fontsize',20)
ylabel('Probability','fontsize',20)               
set(gca,'fontsize',20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% METHOD 1: averaging %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized distance error
figure(3)
normBins = linspace(min(allErr),max(allErr),100);
normDistSubj = NaN(max(subj),100);
for s = 1:max(subj)  
    normDistSubj(s,:) = [histc(allErr(subj==s),normBins)/sum(subj==s)]';    
end
bar(normBins,mean(normDistSubj),1,'facecolor',[.5 .5 .5])
set(gca,'xlim',[0 1])
grid on
xlabel('Norm. Distance Error','fontsize',20)
ylabel('Probability','fontsize',20)               
set(gca,'fontsize',20)

% euclidean distance error
figure(4)
eucBins = linspace(min(allEuc),max(allEuc),100);
eucDistSubj = NaN(max(subj),100);
for s = 1:max(subj)  
    eucDistSubj(s,:) = [histc(allEuc(subj==s),eucBins)/sum(subj==s)]';    
end
bar(eucBins,mean(eucDistSubj),1,'facecolor',[.5 .5 .5])
set(gca,'xlim',[0 max(allEuc)])
grid on
xlabel('Euc. Distance Error','fontsize',20)
ylabel('Probability','fontsize',20)               
set(gca,'fontsize',20)