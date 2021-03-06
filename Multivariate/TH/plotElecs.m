function plotElecs(subj,duration,offset)

% load electrode localizaiton
tal = getBipolarSubjElecs(subj,1)

if mean(strcmp({tal.eType},'D')) == 1
    return
end


% load events
events = get_sub_events('RAM_TH1',subj);
for k=1:length(events)
    events(k).eegfile=regexprep(events(k).eegfile,'eeg.reref','eeg.noreref');
end

% filter to enc
eventsToUse = strcmp({events.type},'CHEST') & [events.confidence]>=0;
encEvents = events(eventsToUse);

evByPowAll  = NaN(sum(eventsToUse),length(tal));
sessions    = [encEvents.session];
uniqSess    = unique(sessions);
for elec = 1:length(tal)
    fprintf('%s: elec %d of %d.\n',subj,elec,length(tal))
    evByPow  = NaN(length(encEvents),1);
    
    for sess = 1:length(uniqSess)
        
        % compute for this session
        sessInds    = sessions==uniqSess(sess);
        sessEv      = encEvents(sessInds);
        [~,amp]     = gethilbertphase_bipol(tal(elec).channel,sessEv,duration,offset,2000,[4 10],60);
        pow         = amp.^2;
        
        [~,pow] = getphasepow_bipol(tal(15).channel,sessEv,duration,offset,2000,'freqs',freqs,'width',5','filtfreq',[58 62],'filttype','stop'); 
        pow(pow<=0) = eps;
        pow         = log10(pow);
        
        % average across time and zscore
        pow = nanmean(pow,3);
        pow = log10(pow);
        pow = zscore(pow);
        %pow         = zscore(nanmean(pow,3));
        %pow = nanmean(pow,3);
        evByPow(sessInds,:) = pow;
    end
    evByPowAll(:,elec) = evByPow;
end



zpow = [];
for sess = 1:length(uniqSess)
    sessInds    = sessions==uniqSess(sess);
        sessEv      = encEvents(sessInds);
    [~,sess_pow] = getphasepow_bipol(tal(15).channel,sessEv,duration,offset,2000,'freqs',freqs,'width',5','filtfreq',[58 62],'filttype','stop'); 
    MU_POW = nanmean(nanmean(sess_pow,3));
    STD_POW = nanstd(nanmean(sess_pow,3));
    [n_ev,n_freq,n_time] = size(sess_pow);
    pat = [];
    for f = 1:n_freq
        all_mu = repmat(MU_POW(f),n_ev,n_time);
        all_std = repmat(STD_POW(f),n_ev,n_time);
        pat(:,f,:) = (squeeze(sess_pow(:,f,:))-all_mu)./all_std;
    end
    zpow = cat(1,zpow,pat);
end

% recalled
recalled = [encEvents.recalled]==1;

recPow = evByPowAll(recalled,:);
nonRecPow = evByPowAll(~recalled,:);
[h,p,c,s] = ttest2(recPow,nonRecPow);

tal_indiv = [tal.indivSurf];
tal_avg = [tal.avgSurf];
if ~isnan(tal_indiv(1).x)
    thisTal = tal_indiv;
else
    thisTal = tal_avg;
end
keyboard

try 
    XYZ = get_xyz(thisTal,'Dykstra');
    if all(all(isnan(XYZ)))
        xyz_field = 'snap';
        XYZ = get_xyz(thisTal, xyz_field);
    end
catch e
    xyz_field = 'snap';
    XYZ = get_xyz(thisTal, 'dural');
end

sigRetain = p' <.05;
posInd = s.tstat' >0;

grpNames = {tal.grpName};
lGrpNames = cellfun(@(x)x(1)=='L',grpNames)';
rGrpNames = cellfun(@(x)x(1)=='R',grpNames)';

l_ind = lGrpNames | (~lGrpNames & ~rGrpNames & XYZ(:,1)<=0);
isDepth = strcmp({tal.eType},'D')';


%% PLOTTING BELOW
%% Left Brain
%l_fig = swagFig; figure(l_fig); ax.L = gca;
clf
surfL = thisTal.path2surfL;
[vL,fL] = read_surf_wrapper(surfL);

plotsurf_wrapper(vL,fL,[.7 .7 .7],0);axis('off');view([-90 0]);
hold on
h = camlight
plot3_wrapper(XYZ(l_ind & ~isDepth,:),55, 'k');

plot3_wrapper(XYZ(l_ind&sigRetain&posInd&~isDepth,:),55, 'r');

plot3_wrapper(XYZ(l_ind&sigRetain&~posInd&~isDepth,:),55, 'b');
set(gca,'visible','off');

%% Right Brain
%r_fig = swagFig; figure(r_fig); ax.R = gca;
surfR = thisTal.path2surfR;
[vR,fR] = read_surf_wrapper(surfR);

plotsurf_wrapper(vR,fR,[.7 .7 .7],0);axis('off');view([90 0]);

hold all
plot3_wrapper(XYZ(~l_ind & ~isDepth,:),55, 'k');
plot3_wrapper(XYZ(~l_ind&sigRetain&posInd& ~isDepth,:),55, 'r');
plot3_wrapper(XYZ(~l_ind&sigRetain&~posInd & ~isDepth,:),55, 'b');
%set(ax.R,'visible','off');
keyboard

