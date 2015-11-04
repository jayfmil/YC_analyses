function YC1_makeCorrPlots(subjs,params)
% function YC1_makeUnivarSubjectReports(subjs,params)
%
% Inputs:
%
%      subjs - cell array of subject strings (default get_subs('RAM_YC1'))
%     params - params structure (default is returned by univarParams)
%  overwrite - boolean. if true, overwrite existing figures
%
%
%
%
%

% if not given, use default params
if ~exist('params','var') || isempty(params)
    params = univarParams();
end


% save directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
region = params.region;
if isempty(region);region = 'all';end
dataDir = fullfile(params.basePath,f(params.useCorrectedPower,y),region);
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

% This loop loads the data from each subject and concatenates all subjects
% in to one large structure
subjDataAll =  [];
for s = 1:length(subjs)                
    
    subj = subjs{s};
    subjFile = fullfile(dataDir,[subj '.mat']);
    if ~exist(subjFile,'file')
        fprintf('Subject data not found for %s.\n',subj)
        continue
    end        

    % load subject data
    subjData = load(subjFile);
    
    % if we haven't done it yet, initialze structure to hold concatenated
    % data from all subjects
    if isempty(subjDataAll)
       fields = fieldnames(subjData.res);
       for f = fields'           
          subjDataAll.(f{1}) = []; 
          subfields = fieldnames(subjData.res.(f{1}));
          for subf = subfields'
              subjDataAll.(f{1}).(subf{1}) = [];
          end
       end
    end    
    
    % merge current subject data into larger struct
    for f = fields' 
        subjDataAll.(f{1}) = mergestruct(subjDataAll.(f{1}),subjData.res.(f{1}));
    end    
end

% now that all the data is loaded, plot stuff..
% bar plot with corr coeff for each freq range
meanR = NaN(1,length(fields));
semR  = NaN(1,length(fields));
pR    = NaN(1,length(fields));

meanRinner = NaN(1,length(fields));
semRinner  = NaN(1,length(fields));

meanRouter = NaN(1,length(fields));
semRouter  = NaN(1,length(fields));

for f = 1:length(fields)
    % note: I'm flipping the sign so that positive correlations indicate
    % better performance
    meanR(f)  = nanmean(-subjDataAll.(fields{f}).r);
    [~,pR(f)] = ttest(-subjDataAll.(fields{f}).r);
    semR(f)   = nanstd(subjDataAll.(fields{f}).r)/sqrt(length(subjDataAll.(fields{f}).r)-1);
    
    % center vertical region
    meanRinner(f)  = nanmean(-subjDataAll.(fields{f}).rInner);
    semRinner(f)   = nanstd(subjDataAll.(fields{f}).rInner)/sqrt(length(subjDataAll.(fields{f}).rInner)-1);

    % outer vertical regions
    meanRouter(f)  = nanmean(-subjDataAll.(fields{f}).rOuter);
    semRouter(f)   = nanstd(subjDataAll.(fields{f}).rOuter)/sqrt(length(subjDataAll.(fields{f}).rOuter)-1);
end

% Figure 1: average correlation for each freq band for all trials
freqs=num2cell(params.freqBins,2);
fun = @(x) strcat(num2str(x(1)),'-',num2str(x(2)));
freqStr = cellfun(fun,freqs,'uniformoutput',false);

figure(1)
clf
bar(find(meanR>0),meanR(meanR>0),'linewidth',4,'FaceColor',[150 23 31]/255);
hold on
bar(find(meanR<0),meanR(meanR<0),'linewidth',4,'FaceColor',[61 89 171]/255);
errorbar(1:f,meanR,semR*1.96,'k','linewidth',4','linestyle','none')
set(gca,'xtick',1:f)
set(gca,'xticklabel',fields')
grid on
set(gca,'gridlinestyle',':')
ylabel('Mean Pearson Coef.','fontsize',26);
set(gca,'fontsize',26)
set(gca,'xlim',[0 length(meanR)+1])
% set(gca,'xticklabel',{'1-3','3-9','40-70','70-200'})
set(gca,'xticklabel',freqStr)
fname = fullfile(figDir,'corrByFreq.eps');
print('-depsc2','-loose',fname)

% Figure 2: average correlation for each freq band for center trials
figure(2)
clf
bar(find(meanRinner>0),meanRinner(meanRinner>0),'linewidth',4,'FaceColor',[150 23 31]/255);
hold on
bar(find(meanRinner<0),meanRinner(meanRinner<0),'linewidth',4,'FaceColor',[61 89 171]/255);
errorbar(1:f,meanRinner,semRinner*1.96,'k','linewidth',4','linestyle','none')
set(gca,'xtick',1:f)
set(gca,'xticklabel',freqStr)
grid on
set(gca,'gridlinestyle',':')
ylabel('Mean Pearson Coef.','fontsize',26);
set(gca,'fontsize',26)
set(gca,'xlim',[0 length(meanR)+1])
fname = fullfile(figDir,'corrByFreqInner.eps');
print('-depsc2','-loose',fname)

% Figure 3: average correlation for each freq band for outer trials
figure(3)
clf
bar(find(meanRouter>0),meanRouter(meanRouter>0),'linewidth',3,'FaceColor',[150 23 31]/255);
hold on
bar(find(meanRouter<0),meanRouter(meanRouter<0),'linewidth',3,'FaceColor',[61 89 171]/255);
errorbar(1:f,meanRouter,semRouter*1.96,'k','linewidth',3','linestyle','none')
set(gca,'xtick',1:f)
set(gca,'xticklabel',freqStr)
grid on
set(gca,'gridlinestyle',':')
ylabel('Mean Pearson Coef.','fontsize',26);
set(gca,'fontsize',26)
set(gca,'xlim',[0 length(meanR)+1])
fname = fullfile(figDir,'corrByFreqOuter.eps');
print('-depsc2','-loose',fname)
keyboard

function sout = mergestruct(struct1,struct2)

if isempty(struct1) & ~isempty(struct2)
    sout = struct2;
    return
elseif ~isempty(struct1) & isempty(struct2)
    sout = struct1;
    return
end

fields1 = fieldnames(struct1);
fields2 = fieldnames(struct2);

if isequal(fields1,fields2)
    sout = cell2struct(fields1,fields1,1);
    for f = 1:length(fields1)
        if ndims(struct1.(fields1{f})) > 2 | ndims(struct2.(fields1{f})) > 2
            sout.(fields1{f}) = cat(3,struct1.(fields1{f}),struct2.(fields1{f}));       
        else
            sout.(fields1{f}) = horzcat(struct1.(fields1{f}),struct2.(fields1{f}));       
        end
    end
end


