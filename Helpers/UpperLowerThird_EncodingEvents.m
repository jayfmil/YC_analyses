function eventMask = UpperLowerThird_EncodingEvents(events)

perfErr = [events.testError];
lowThresh = prctile(perfErr,33);
highThresh = prctile(perfErr,67);

% filter to just encoding events
eventMask = (strcmp({events.type},'NAV_LEARN') | strcmp({events.type},'NAV_LEARN_HARD')) & ((perfErr <= lowThresh) | (perfErr >= highThresh));

% determine thresholds for upper third of performance and lower third
% perfErr = [events(allEncoding).testError];
% lowThresh = prctile(perfErr,33);
% highThresh = prctile(perfErr,67);
% 
% % restrict to just upper or lower third events
% eventMask = (perfErr <= lowThresh) | (perfErr >= highThresh);