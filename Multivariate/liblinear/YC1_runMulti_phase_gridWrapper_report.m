function YC1_runMulti_phase_gridWrapper_report(subjs,params)

basePath = '/scratch/jfm2/YC1/multi/watrous_gridWrapper';
basePath = '/scratch/jfm2/YC1/multi/pow8_gridWrapper';
basePath = '/scratch/jfm2/YC1/multi/pow8_gridWrapper_fixed';
% basePath = '/scratch/jfm2/YC1/multi/pow8_gridWrapper_fixed_block'

figDir = fullfile(basePath,'figs');
if ~exist(figDir,'dir')
    mkdir(figDir);
end

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end

% use default params if none given
if ~exist('params','var') || isempty(params)
    params = multiParams();
end
params.crossValStrictness = 1;
params.useKfold = 0;
params.cvField  = 'session';



aucs_all = NaN(4,2,2,3,3,length(subjs));

% loop over subjets
for s = 1:length(subjs)
    
    % loop over time bin settings
    timeBins = {[1 20;21 40;41 60;61 80;81 100],[1 20],[81 100],[1 100]};
    for t = 1:length(timeBins)
        params.timeBins = timeBins{t};
        
        % loop over regularization method
        normTypes = {'L1','L2'};
        for l = 1:2
            params.normType = normTypes{l};
            
            % loop over inner cross val scheme
            nestedCvFields = {'session','blocknum'};
            for cv = 1:1
                params.nestedCvField = nestedCvFields{cv};
                
                % loop over encoding period selection
                encPeriods = {'first','second','both'};
                for enc = 1:3
                    params.encPeriod = encPeriods{enc};
                    
                    % make save directory
                    pathExt = sprintf('t-%d_l-%d_cv-%d_enc-%d',t,l,cv,enc);
                    saveDir = fullfile(basePath,pathExt);
                    
                    % lastly, loop over phase, power, or powre and phase
                    % together (0=power,1=phase,2=both)
                    usePhase = [0,1,2];
                    for p = 1:3
                        params.usePhase = usePhase(p);
                        fname = fullfile(saveDir,[subjs{s} '_lasso_pow.mat']);
                        if params.usePhase==1
                            fname = fullfile(saveDir,[subjs{s} '_lasso_phase.mat']);
                        elseif params.usePhase==2
                            fname = fullfile(saveDir,[subjs{s} '_lasso_powphase.mat']);
                        end
                        if ~exist(fname,'file')                                                        
                            continue
                        else
                            load(fname)
                            aucs_all(t,l,cv,enc,p,s) = AUC;
                        end
                    end
                end
            end
        end
    end
end

sessSubjs   = ~isnan(squeeze(aucs_all(1,1,1,1,1,:)));
blockSubjs  = ~isnan(squeeze(aucs_all(1,1,2,1,1,:)));



aucsSessL1_phase = nanmean(squeeze(aucs_all(:,1,1,:,2,:)),3);
aucsSessL2_phase = nanmean(squeeze(aucs_all(:,2,1,:,2,:)),3);
aucsBlockL1_phase = nanmean(squeeze(aucs_all(:,1,2,:,2,:)),3);
aucsBlockL2_phase = nanmean(squeeze(aucs_all(:,2,2,:,2,:)),3);

aucsSessL1_pow = nanmean(squeeze(aucs_all(:,1,1,:,1,:)),3);
aucsSessL2_pow = nanmean(squeeze(aucs_all(:,2,1,:,1,:)),3);
aucsBlockL1_pow = nanmean(squeeze(aucs_all(:,1,2,:,1,:)),3);
aucsBlockL2_pow = nanmean(squeeze(aucs_all(:,2,2,:,1,:)),3);
aucsBlockL1_matched_pow = nanmean(squeeze(aucs_all(:,1,2,:,1,sessSubjs)),3);
aucsBlockL2_matched_pow = nanmean(squeeze(aucs_all(:,2,2,:,1,sessSubjs)),3);


aucsSessL1_both = nanmean(squeeze(aucs_all(:,1,1,:,3,:)),3);
aucsSessL2_both = nanmean(squeeze(aucs_all(:,2,1,:,3,:)),3);
aucsBlockL1_both = nanmean(squeeze(aucs_all(:,1,2,:,3,:)),3);
aucsBlockL2_both = nanmean(squeeze(aucs_all(:,2,2,:,3,:)),3);
keyboard

data = {aucsSessL1_pow,aucsSessL2_pow,aucsBlockL1_pow,aucsBlockL2_pow,aucsBlockL1_matched_pow,aucsBlockL2_matched_pow};
titleNormtypes = {'L1','L2','L1','L2','L1','L2'};
titleCV = {'Session','Session','Trial','Trial','Trial','Trial'};
titleCounts = [sum(sessSubjs) sum(sessSubjs) sum(blockSubjs) sum(blockSubjs) sum(sessSubjs) sum(sessSubjs)];
fnames = {'aucsSessL1_pow','aucsSessL2_pow','aucsBlockL1_pow','aucsBlockL2_pow','aucsBlockL1_matched_pow','aucsBlockL2_matched_pow'};

figure(1)
for i = 1:6
    clf
    plotData = data{i};
    imagesc(plotData)
    set(gca,'clim',[.45 .6])
    set(gca,'xtick',1:3)
    set(gca,'xticklabel',{'Enc 1','Enc 2','Both'})
    xlabel('Encoding Period','fontsize',20)
    set(gca,'ytick',1:4)
    set(gca,'yticklabel',{'1s Bins','Start','End','5s Avg'})
    ylabel('Time Period','fontsize',20)
    set(gca,'fontsize',20)
    colorbar
    titleStr = sprintf('%s %s CV, %d subjects, max AUC: %.2f',titleNormtypes{i},titleCV{i},titleCounts(i),max(plotData(:)));
    title(titleStr);
    set(gca,'titlefontweight','normal');
    fname = fullfile(figDir,fnames{i});
    print('-depsc2','-loose',fname)    
end























