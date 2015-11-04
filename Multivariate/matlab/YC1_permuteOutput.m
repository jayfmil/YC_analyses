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

if params.modelEachTime
    nTimes = size(params.timeBins,1);
else
    nTimes = 1;
end

AUC  = NaN(1,nTimes);
perf = NaN(1,nTimes);
r    = NaN(1,nTimes); 
% loop over each time bin

if params.modelEachTime
    nTimes = size(params.timeBins,1);
else
    nTimes = 1;
end

for t = 1:nTimes
    
    % get the actual predicted values
    yPred     = vertcat(subjData.res(t).yPred{:});
    
    if isfield(params,'encPeriod') && (strcmpi(params.encPeriod,'combined')  || strcmpi(params.encPeriod,'average') || strcmpi(params.encPeriod,'first') || strcmpi(params.encPeriod,'second')) 
        randOrder = randperm(length(yPred));
    else
        randOrder = randperm(length(yPred)/2);
    end
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