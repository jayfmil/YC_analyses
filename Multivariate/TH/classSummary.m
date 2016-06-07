function classSummary

baseDir      = '/scratch/jfm2/TH1/multi';
classSchemes = {'gridsearch_varWindowSize','gridsearch_varWindowSize_conf','gridsearch_varWindowSize_medianThresh','gridsearch_varWindowSize_lowWrong_median'};
dataBin      = 'w_9_b_16_f_3_r_1';
fullDirs     = fullfile(baseDir,classSchemes,dataBin);



subjs = get_subs('RAM_TH1');
aucs = NaN(length(subjs),length(fullDirs));
ps   = NaN(length(subjs),length(fullDirs));
for i = 1:length(fullDirs)
    load(fullfile(fullDirs{i},[subjs{1} '_class_pow.mat']));
    [aucs(:,i),ps(:,i)] = TH_plotClassRes(subjs,params,0,1);
end
keyboard

is_loso = false(length(subjs),1);
for s = 1:length(subjs)
    events = get_sub_events('RAM_TH1',subjs{s});
    if length(unique([events.session])) > 1
        is_loso(s) = 1;
    end
end