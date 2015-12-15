function YC1_makeElecTable(subjs)


% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

regions = {'Hippform','Hipp','CA1','CA2','CA3','DG','Sub','EC'};
table = NaN(length(subjs),length(regions));

for s = 1:length(subjs)
    
    try
        tal = getBipolarSubjElecs(subjs{s});
        
        for r = 1:length(regions)
            
            talRegion = filterTalByRegion(tal,regions{r});
            table(s,r) = length(talRegion);
        end        
    end
end
keyboard


