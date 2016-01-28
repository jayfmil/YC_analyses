function figs = YC1_plotPowerSpect(subjs,params)
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
y2 = {'','average'};
region = params.region;
if isempty(region);region = 'all';end
dataDir = fullfile(params.basePath,f(params.useCorrectedPower,y),region,f(params.averageRegion,y2));
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

% structre to hold figure paths
figs = [];

for s = 1:length(subjs)                
    
    subj = subjs{s};
    subjFile = fullfile(dataDir,[subj '.mat']);    
    if ~exist(subjFile,'file')
        fprintf('Subject data not found for %s.\n',subj)
        continue
    end        
    subjData = load(subjFile);
    figsSubj = [];
    
    inner = subjData.isInnerBand;    
    for elec = 1:size(subjData.fullPowerSpect,2);
       
        if length(subjData.tal(elec).channel) == 1
            chanStr = num2str(subjData.tal(elec).channel);
        else
            chanStr = [num2str(subjData.tal(elec).channel(1)) '-' num2str(subjData.tal(elec).channel(2))];
        end
        
        % bin by inner or outer area
        elecData = squeeze(subjData.fullPowerSpect(:,elec,:));
        innerData = elecData(inner,:);
        outerData = elecData(~inner,:);
        
        % plot power spectrum
        clf
        handle = plotLogLogPowSpect(innerData,outerData,subjData.freqs,24,20);        
        hold on
        [h,p,c,st] = ttest2(innerData,outerData);
        pthresh = .05/length(subjData.freqs);
        pthresh = .05;
        if any(p < pthresh)
            y = ylim;
            set(gca,'ylim',[y(1) y(2)+.5])
            y = max(ylim)-.25;
            colors = {[226,55,67]/255,[50,124,203]/255};
            sig = p < pthresh;
            if any(st.tstat > 0 & sig)
                plot(log10(subjData.freqs(st.tstat > 0 & sig)),y,'.','markersize',20,'color',colors{1});
            end
            if any(st.tstat < 0 & sig)
                plot(log10(subjData.freqs(st.tstat < 0 & sig)),y,'.','markersize',20,'color',colors{2});
            end            
        end        
        legend(handle,'Inner','Outer')
        titleStr = sprintf('%s: %s - %s',subjs{s},chanStr,subjData.locTags{elec});
        title(titleStr);
        set(gca,'TitleFontWeight','normal')
        
        % save figure
        fname = fullfile(figDir,sprintf('%s_pSpect_%s.eps',subjs{s},chanStr));
        print('-depsc2','-loose',fname);
        figsSubj(elec).channel = chanStr;       
        figsSubj(elec).loc = subjData.locTags{elec};
        figsSubj(elec).fname = fname;                
    end    
    figs.(subjs{s}) = figsSubj;
end










