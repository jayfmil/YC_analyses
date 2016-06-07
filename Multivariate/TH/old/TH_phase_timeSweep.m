function TH_phase_timeSweep(subjs)


% p = TH_multiParams
% % params.cvField = 'session';
% params.timeBinLabels = {''};
% params.powerPath = '/scratch/jfm2/power8freqs/';
% basePath = '/scratch/jfm2/TH1/multi/acrossTrial';
% params.normType = 'L2';
% params.freqBins = [];

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_TH1');
end


Cs = logspace(log10(1e-6),log10(1e4),22);
Ts = 25:125;
% aucs = NaN(length(subjs),length(Ts));
% p = NaN(length(subjs),length(Ts));
% for s = 1:length(subjs)
%     
%     
%     fprintf('Processing %s.\n',subjs{s})
%     aucs_subj = NaN(1,length(Ts));
%     parfor t = 1:length(Ts);
%         params = TH_multiParams();
%         params.Cs = params.Cs(7);
%         %         params.Cs = .1389;
%         params.powerPath = '/scratch/jfm2/power4bins_hilbert/';
%         params.basePath = '/scratch/jfm2/TH1/multi/acrossTrial_timeSweep';
%         params.normType = 'L2';
%         params.freqBins = [];
%         params.saveOutput = 0;
%         params.timeBins = Ts(t);
%         params.freqBins = [];
%         params.usePhase = 1;
%         
%         saveDir = fullfile(params.basePath,['_',num2str(t)]);
%         if ~exist(saveDir,'dir')
%             mkdir(saveDir)
%         end
%         %         clear AUC res
%         %         if exist(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']),'file')
%         %             load(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']));
%         %         else
%         AUC = TH1_refactor_phase(subjs{s},params,saveDir);
%         %         end
%         %         if ~isempty(AUC)
%         aucs_subj(t) = AUC;
%         %         end
%     end
%     aucs(s,:) = aucs_subj;
%     
% end
keyboard
for s = 1:length(subjs)
    
    
    fprintf('Processing %s.\n',subjs{s})    
    for t = 1:length(Ts);
        
        aucs_t = NaN(1,100);
        parfor iter = 1:100
            params = TH_multiParams();
            params.Cs = params.Cs(7);
            %         params.Cs = .1389;
            params.powerPath = '/scratch/jfm2/power4bins_hilbert/';
            params.basePath = '/scratch/jfm2/TH1/multi/acrossTrial_timeSweep';
            params.normType = 'L2';
            params.freqBins = [];
            params.saveOutput = 0;
            params.timeBins = Ts(t);
            params.freqBins = [];
            params.usePhase = 1;
            params.loadPower = 1;
            params.doPermute = 1;
            
            saveDir = fullfile(params.basePath,['_',num2str(t)]);
            %         clear AUC res
            %         if exist(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']),'file')
            %             load(fullfile(params.basePath,[subjs{s} '_lasso_pow.mat']));
            %         else
            aucs_t(iter) = TH1_refactor_phase(subjs{s},params,saveDir);
            %         end
            %         if ~isempty(AUC)
            %         end
        end
        p(s,t) = mean(aucs(s,t) < aucs_t)
    end    
    
end

keyboard
basePath = '/scratch/jfm2/TH1/multi/acrossTrial_timeSweep';
fname = fullfile(basePath,'aucs_timeSweep.mat');
save(fname,'aucs','Ts','subjs','p');
keyboard











