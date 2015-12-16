function events = addExtraYCFields(events)
% add testError field
% add inner field (1 = inner region, 0 = outer region)

testInd = strcmp({events.type},'NAV_TEST');
recEvents = events(testInd);
[events.testError] = deal(NaN);
[events.recalled] = deal(NaN);
[events.inner] = deal(NaN);
[events.sessionHalf] = deal(0);
[events.withinItemCount] = deal(0);
sessVec = [events.session];
trialVec = [events.blocknum];

uniqSess = unique(sessVec);
sessPerc = NaN(1,length(sessVec));

for sess = 1:length(uniqSess)
    sessInd = sessVec == uniqSess(sess);
    sessPerc(sessInd) = [1:sum(sessInd)]./sum(sessInd);
    
    items = unique(trialVec(sessInd));
    for item = 1:length(items)
        ind = trialVec == items(item) & sessInd;
        
        for c = 1:sum(ind)
            indTmp = find(ind);
            events(indTmp(c)).withinItemCount = c;
        end
    end
    
end
[events(sessPerc>.5).sessionHalf] = deal(1);


for rec = 1:length(recEvents);
    session = recEvents(rec).session;
    trial = recEvents(rec).blocknum;
    err = recEvents(rec).respPerformanceFactor;
    ind = sessVec == session & trialVec == trial;
    [events(ind).testError] = deal(err);
    [events(ind).inner] = deal(abs(recEvents(rec).objLocs(1)) < 568/30 && abs(recEvents(rec).objLocs(2)) < 7);
end



