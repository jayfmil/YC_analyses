function plot_R1054_perm



baseDir = '/scratch/jfm2/YC1/multi/R1054J_new';
randNeural = load(fullfile(baseDir,'randX','randDist.mat'));
randResp = load(fullfile(baseDir,'permute','randDist.mat'));
randNeural_samp = load(fullfile(baseDir,'randX_sample','randDist.mat'));
randResp_train = load(fullfile(baseDir,'permute_train','randDist.mat'));

figure
clf
[n,x] = hist(randNeural.perfs,25);
bar(x,n/sum(n),'hist','w');
hold on
y = binopdf(1:48,48,.5);
plot([1:48]/48,y)
grid on
set(gca,'gridlinestyle',':')
xlabel('Prob. Correct','fontsize',20);
ylabel('Probability','fontsize',20);
set(gca,'fontsize',20)
print('-depsc2','-loose','~/R1054J_randNeural')


figure
clf
[n,x] = hist(randResp.perfs,25);
bar(x,n/sum(n),'hist','w');
hold on
y = binopdf(1:48,48,.5);
plot([1:48]/48,y)
grid on
set(gca,'gridlinestyle',':')
xlabel('Prob. Correct','fontsize',20);
ylabel('Probability','fontsize',20);
set(gca,'fontsize',20)
print('-depsc2','-loose','~/R1054J_randResp_sample')

figure
clf
[n,x] = hist(randNeural_samp.perfs,25);
bar(x,n/sum(n),'hist','w');
hold on
y = binopdf(1:48,48,.5);
plot([1:48]/48,y)
grid on
set(gca,'gridlinestyle',':')
xlabel('Prob. Correct','fontsize',20);
ylabel('Probability','fontsize',20);
set(gca,'fontsize',20)
print('-depsc2','-loose','~/R1054J_randNeural_samp')


figure
clf
[n,x] = hist(randResp_train.perfs,25);
bar(x,n/sum(n),'hist','w');
hold on
y = binopdf(1:48,48,.5);
plot([1:48]/48,y)
grid on
set(gca,'gridlinestyle',':')
xlabel('Prob. Correct','fontsize',20);
ylabel('Probability','fontsize',20);
set(gca,'fontsize',20)
print('-depsc2','-loose','~/R1054J_randResp_train')