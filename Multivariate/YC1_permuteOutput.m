function [perf,AUC,r] = YC1_permuteOutput(subj,params,saveDir)
%function [perf,auc,r] = YC1_permuteOutput(subj,params,saveDir)
%
%
%
perf    = [];
AUC     = [];
r       = [];

% load subject data
subjData = load(fullfile(saveDir,[subj '_lasso.mat']));

AUC  = NaN(1,size(params.timeBins,1));
perf = NaN(1,size(params.timeBins,1));
r    = NaN(1,size(params.timeBins,1));
% loop over each time bin
for t = 1:size(params.timeBins,1)
    
    % get the actual predicted values
    yPred     = vertcat(subjData.res(t).yPred{:});
    
    randOrder = randperm(length(subjData.res(t).yPred));
    if ~isfield(params,'encPeriod') || strcmpi(params.encPeriod,'both')        
        randOrder = [randOrder;randOrder];
        randOrder = randOrder(:);
    end
        
    % permute the responses
    Y = subjData.Y;
    yPerm = Y(randOrder);
    
    % calculate new perf
    perf(t) = mean((yPerm >= mean(Y)) == Y);
    
    % calculate new AUC
    [~,~,~,AUC(t)] = perfcurve(yPerm,yPred,true);
    
end