function YC1_perfByCs(subjs,basePath)



numTs   = 4;
numLs   = 2;
numCvs  = 2;
numEncs = 3;
numCs   = 22;

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC1');
end


aucs_all = NaN(numTs,numLs,numCvs,numEncs,numCs,length(subjs));

% loop over subjets
for s = 1:length(subjs)
    s
    % loop over time bins
    for t = 1:numTs
        
        % loop over regularization method
        for l = 1:numLs
            
            % loop over inner cross val scheme
            for cv = 1:numCvs
                
                % loop over encoding period selection
                for enc = 1:numEncs
                    
                    % get save directory
                    pathExt = sprintf('t-%d_l-%d_cv-%d_enc-%d',t,l,cv,enc);
                    
                    for c = 1:numCs;
                        
                        dataDir = fullfile(basePath,pathExt,sprintf('Cs_%d',c));
                        fname = fullfile(dataDir,[subjs{s} '_lasso_pow.mat']);
                        if ~java.io.File(fname).exists;%~exist(fname,'file')
                            continue
                        end                        
                        load(fname,'AUC')
                        aucs_all(t,l,cv,enc,c,s) = AUC;                                                
                    end                    
                end
            end
        end
    end
end



mean_l1_sess  = NaN(numTs,numCs);
mean_l2_sess  = NaN(numTs,numCs);
mean_l1_block = NaN(numTs,numCs);
mean_l2_block = NaN(numTs,numCs);

e_l1_sess  = NaN(numTs,numCs);
e_l2_sess  = NaN(numTs,numCs);
e_l1_block = NaN(numTs,numCs);
e_l2_block = NaN(numTs,numCs);

p_l1_sess  = NaN(numTs,numCs);
p_l2_sess  = NaN(numTs,numCs);
p_l1_block = NaN(numTs,numCs);
p_l2_block = NaN(numTs,numCs);


for t = 1:numTs
    
    % get means of l1/l2 and block/sess for each C
    mean_l1_sess(t,:)  = squeeze(nanmean(aucs_all(t,1,1,3,:,:),6))';
    mean_l2_sess(t,:)  = squeeze(nanmean(aucs_all(t,2,1,3,:,:),6))';
    mean_l1_block(t,:) = squeeze(nanmean(aucs_all(t,1,2,3,:,:),6))';
    mean_l2_block(t,:) = squeeze(nanmean(aucs_all(t,2,2,3,:,:),6))';
    
    % get standard errors
    e_l1_sess(t,:)  = squeeze(nanstd(aucs_all(t,1,1,3,:,:),[],6))';
    nSubj           = sum(squeeze(~isnan(aucs_all(t,1,1,3,:,:)))',1);
    e_l1_sess(t,:) = e_l1_sess(t,:)./sqrt(nSubj);
    
    e_l2_sess(t,:)  = squeeze(nanstd(aucs_all(t,2,1,3,:,:),[],6))';
    nSubj           = sum(squeeze(~isnan(aucs_all(t,2,1,3,:,:)))',1);
    e_l2_sess(t,:) = e_l2_sess(t,:)./sqrt(nSubj);
    
    e_l1_block(t,:) = squeeze(nanstd(aucs_all(t,1,2,3,:,:),[],6))';
    nSubj           = sum(squeeze(~isnan(aucs_all(t,1,2,3,:,:)))',1);
    e_l1_block(t,:) = e_l1_block(t,:)./sqrt(nSubj);
    
    e_l2_block(t,:) = squeeze(nanstd(aucs_all(t,2,2,3,:,:),[],6))';    
    nSubj           = sum(squeeze(~isnan(aucs_all(t,2,2,3,:,:)))',1);
    e_l2_block(t,:) = e_l2_block(t,:)./sqrt(nSubj);
    
    
    % run ttests vs .5
    [h,p,c,s] = ttest(squeeze(aucs_all(t,1,1,3,:,:))',.5);
    p_l1_sess(t,:) = p;
    
    [h,p,c,s] = ttest(squeeze(aucs_all(t,2,1,3,:,:))',.5);
    p_l2_sess(t,:) = p;
    
    [h,p,c,s] = ttest(squeeze(aucs_all(t,1,2,3,:,:))',.5);
    p_l1_block(t,:) = p;
    
    [h,p,c,s] = ttest(squeeze(aucs_all(t,2,2,3,:,:))',.5);
    p_l2_block(t,:) = p;    
    
    % plot
    
end
keyboard
figDir = fullfile(basePath,'figs');
if ~exist(figDir,'dir')
    mkdir(figDir);
end

cL1 = logspace(log10(1e-2),log10(1e4),22);
cL2 = logspace(log10(1e-6),log10(1e4),22);

close all
for t = 1:4        
    
    %%% -------------------------------------------------------------------
    % Sess
    
    % L1 Sess
    figure(1)
    subplot(4,1,t)
    errorbar(1:numCs,mean_l1_sess(t,:),e_l1_sess(t,:)*1.96,'-k','linewidth',3);
    set(gca,'ylim',[.45 .6])
    set(gca,'xtick',[1 10 22])
    grid on
    hold on
    plot(xlim,[.5 .5],'--k','linewidth',2)
    
    sig = find(p_l1_sess(t,:)<.05);
    for s = 1:length(sig)                
        errorbar(sig(s),mean_l1_sess(t,sig(s)),e_l1_sess(t,sig(s))*1.96,'r','linewidth',3);
    end
    if t == 4
        set(gca,'xticklabel',cL1([1 10 22]))    
        set(gca,'gridlinestyle',':');
        ylabel('AUC','fontsize',20)
        xlabel('Cost Parameter','fontsize',20);
    else
        set(gca,'xticklabel','')
    end
    set(gca,'fontsize',20)  
    print('-depsc2','-loose',fullfile(figDir,'L1_sess'))
    
    % L2 Sess
    figure(2)
    subplot(4,1,t)
    errorbar(1:numCs,mean_l2_sess(t,:),e_l2_sess(t,:)*1.96,'-k','linewidth',3);
    set(gca,'ylim',[.45 .6])
    set(gca,'xtick',[1 10 22])
    grid on
    hold on
    plot(xlim,[.5 .5],'--k','linewidth',2)
    
    sig = find(p_l2_sess(t,:)<.05);
    for s = 1:length(sig)                
        errorbar(sig(s),mean_l2_sess(t,sig(s)),e_l2_sess(t,sig(s))*1.96,'r','linewidth',3);
    end
    if t == 4
        set(gca,'xticklabel',cL2([1 10 22]))    
        set(gca,'gridlinestyle',':');
        ylabel('AUC','fontsize',20)
        xlabel('Cost Parameter','fontsize',20);
    else
        set(gca,'xticklabel','')
    end
    set(gca,'fontsize',20)    
    print('-depsc2','-loose',fullfile(figDir,'L2_sess'))
    %%% -------------------------------------------------------------------
    
    %%% -------------------------------------------------------------------
    % Block    
    % L1 Block
    figure(3)
    subplot(4,1,t)
    errorbar(1:numCs,mean_l1_block(t,:),e_l1_block(t,:)*1.96,'-k','linewidth',3);
    set(gca,'ylim',[.45 .6])
    set(gca,'xtick',[1 10 22])
    grid on
    hold on
    plot(xlim,[.5 .5],'--k','linewidth',2)
    
    sig = find(p_l1_sess(t,:)<.05);
    for s = 1:length(sig)                
        errorbar(sig(s),mean_l1_block(t,sig(s)),e_l1_block(t,sig(s))*1.96,'r','linewidth',3);
    end
    if t == 4
        set(gca,'xticklabel',cL1([1 10 22]))    
        set(gca,'gridlinestyle',':');
        ylabel('AUC','fontsize',20)
        xlabel('Cost Parameter','fontsize',20);
    else
        set(gca,'xticklabel','')
    end
    set(gca,'fontsize',20)  
    print('-depsc2','-loose',fullfile(figDir,'L1_block'))
    
    % L2 Block
    figure(4)
    subplot(4,1,t)
    errorbar(1:numCs,mean_l2_block(t,:),e_l2_block(t,:)*1.96,'-k','linewidth',3);
    set(gca,'ylim',[.45 .6])
    set(gca,'xtick',[1 10 22])
    grid on
    hold on
    plot(xlim,[.5 .5],'--k','linewidth',2)
    
    sig = find(p_l2_sess(t,:)<.05);
    for s = 1:length(sig)                
        errorbar(sig(s),mean_l2_block(t,sig(s)),e_l2_block(t,sig(s))*1.96,'r','linewidth',3);
    end
    if t == 4
        set(gca,'xticklabel',cL2([1 10 22]))    
        set(gca,'gridlinestyle',':');
        ylabel('AUC','fontsize',20)
        xlabel('Cost Parameter','fontsize',20);
    else
        set(gca,'xticklabel','')
    end
    set(gca,'fontsize',20)   
    print('-depsc2','-loose',fullfile(figDir,'L2_block'))
    %%% -------------------------------------------------------------------
    
end


















