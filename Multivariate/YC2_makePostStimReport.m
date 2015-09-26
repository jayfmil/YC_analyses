function YC2_makePostStimReport(subjs,params,overwrite)
% function YC2_makePostStimReport(subjs,params,overwrite)
%
% Inputs:
%
%      subjs - cell array of subject strings (default get_subs('RAM_YC1'))
%     params - params structure (default is returned by multiParams)
%  overwrite - boolean. if true, overwrite existing figures
%


% UPDATE THIS
% Make a report of the classifier performance for each subject in YC1 that
% has classification and chance clasification data. Also make a group
% average report.
%
% For each subject, figures are:
%
%   - Classifier performance over time for % correct and for AUC
%   - Recall change as a function of classifier output for the best
%     performing  time bin
%   - Histogram of the patient's behavioral perfomance
%
%
% For the group report, figures are quartile plots for each time bin,
% accuracy historograms on the subject level for each time bin, and AUC
% histograms on the subject level for each time bin.
%
%  Reports are saved to <location of subject
%  data>/reports/lassoChance_report.pdf and group_lassoChance_report.pdf

% if not given, use default params
if ~exist('params','var') || isempty(params)
    params = multiParams();
end

% overwrite exisiting figures?
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = true;
end

% tex directory
f = @(x,y) y{double(x)+1};
y = {'OrigPower','CorrectedPower'};
dataDir = fullfile(params.basePath,f(params.useCorrectedPower,y),'YC2_postStimFFT_both_best')
saveDir = fullfile(dataDir,'reports');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% get list of YC subjects
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs('RAM_YC2');
end
% subjs = subjs(~strcmp(subjs,'R1061T'))

% figure directory
figDir = fullfile(saveDir,'figs');
if ~exist('figDir','dir')
    mkdir(figDir)
end

% i would prallocate but I don't have many regions each subject has
deltaRR      = [];
deltaRR_bin  = [];
deltaEE      = [];
deltaEE_Prob = [];
region       = {};
subjects     = {};
yc1Score     = [];
numNZweights  = [];

% hold the permutaiton results
deltaRR_perm           = [];
deltaRR_bin_perm       = [];
deltaEE_perm           = [];
deltaEE_Prob_perm      = [];
postProbNonStim_perm   = [];
postProbStim_perm      = [];
logOddsNonStim_perm    = [];
logOddsStim_perm       = [];
stimPerf_perm          = [];
nonStimPerf_perm       = [];

% will hold figure paths for latex report
figs = [];
for s = 1:length(subjs)
    subj  = subjs{s};
    fname = fullfile(dataDir,[subj '_YC2_postStimChange.mat']);
    if ~exist(fname,'file')
        continue
    end
    
    % load subject post stim data
    subjData = load(fname);
    
    nRegions = length(subjData.res);
    for r = 1:nRegions
        deltaRR               = [deltaRR;subjData.res(r).deltaRR];
        deltaRR_bin           = [deltaRR_bin;subjData.res(r).deltaRR_bin];
        deltaEE               = [deltaEE;subjData.res(r).deltaEE];
        numNZweights          = [numNZweights;subjData.res(r).numNZweights];
        deltaEE_Prob          = [deltaEE_Prob;subjData.res(r).deltaEE_Prob];
        deltaRR_perm          = [deltaRR_perm;subjData.res(r).deltaRR_perm];
        deltaRR_bin_perm      = [deltaRR_bin_perm;subjData.res(r).deltaRR_bin_perm];
        deltaEE_perm          = [deltaEE_perm;subjData.res(r).deltaEE_perm];
        deltaEE_Prob_perm     = [deltaEE_Prob_perm;subjData.res(r).deltaEE_Prob_perm]; 
        postProbNonStim_perm  = [postProbNonStim_perm;subjData.res(r).postProbNonStim_perm]; 
        postProbStim_perm     = [postProbStim_perm;subjData.res(r).postProbStim_perm]; 
        logOddsNonStim_perm   = [logOddsNonStim_perm;subjData.res(r).logOddsNonStim_perm]; 
        logOddsStim_perm      = [logOddsStim_perm;subjData.res(r).logOddsStim_perm]; 
        stimPerf_perm         = [stimPerf_perm;subjData.res(r).stimPerf_perm]; 
        nonStimPerf_perm      = [nonStimPerf_perm;subjData.res(r).nonStimPerf_perm]; 
        
        stimAnat = subjData.res(r).stimAnat;
        if isempty(subjData.res(r).stimAnat)
            stimAnat = '';
        end
        region{end+1} = stimAnat;
        subjects{end+1} = subj;
        yc1Score      = [yc1Score;subjData.res(r).yc1Score];                
    end    
end

ec   = ~cellfun('isempty',regexpi(region,['ec']))';
hipp = ~cellfun('isempty',regexpi(region,['ca1|ca2|ca3|dg|sub']))';
mtl  = ~cellfun('isempty',regexpi(region,['amy|phc|prc|BA36|pcg']))';
oth  = ~(ec | hipp | mtl);


deltaEE(abs((deltaEE - nanmean(deltaEE))/nanstd(deltaEE)) > 3) = NaN;
xdata      = {deltaEE_Prob,deltaEE};
ydata      = {deltaRR,deltaRR_bin};
xDataNames  = {'prob','logOdds'};
yDataNames  = {'perfScore','RecChange'};
subjStr    = {'all','good','weights'};
xdata_perm = {deltaEE_Prob_perm,deltaEE_perm};
ydata_perm = {deltaRR_perm,deltaRR_bin_perm};
subjFilter = {@(x) (~isnan(x) & ~isnan(deltaRR)),@(x) (~isnan(x) & ~isnan(deltaRR) & yc1Score>.95),@(x) (~isnan(x) & ~isnan(deltaRR) &~oth)};
% subjFilter = {@(x) (~isnan(x) & ~isnan(deltaRR) &~oth),@(x)
% (~isnan(x) & ~isnan(deltaRR) & yc1Score>.95 &~oth),@(x)
% (~isnan(x) &  numNZweights>5 &~oth)};
subjFilter = {@(x) (~isnan(x) & ~isnan(deltaRR)),@(x) (~isnan(x) & ~isnan(deltaRR) & yc1Score>.95),@(x) (~isnan(x) & ~isnan(deltaRR) &  ~oth)};

%%% Plots are scatter plots of change in classifier encoding estimate post
%%% stim and change in behavioral performance post stim. Four plots:
%
% 1-2: percent change in classifier PROBABILITY for all subjects and just
%      subjects with good YC1 decoding
% 3-4: percent change in classifier log odds for all subjects and just
%      subjects with good YC1 decoding
figure(1)
figs = [];

for xDataType = 1:2
    x = xdata{xDataType};
    for yDataType = 1:2
        y = ydata{yDataType};
        for filterInd = 1:3
            
            fname = sprintf('postStim_%s_%s_%s',xDataNames{xDataType},yDataNames{yDataType},subjStr{filterInd});
            figs.(fname) = fullfile(figDir,fname);
            
            % scatter with different color for different regions
            subjsToUse = subjFilter{filterInd}(x);
            clf
            plot(x(ec&subjsToUse),y(ec&subjsToUse),'.','markersize',40,'color',[0, 0.4470, 0.7410],'linewidth',3)
            hold on
            plot(x(hipp&subjsToUse),y(hipp&subjsToUse),'.','markersize',40,'color',[0.8500, 0.3250, 0.0980],'linewidth',3)
            plot(x(mtl&subjsToUse),y(mtl&subjsToUse),'.','markersize',40,'color',[0.3010, 0.7450, 0.9330],'linewidth',3)
            plot(x(oth&subjsToUse),y(oth&subjsToUse),'.k','markersize',40)
            s = regstats(y(subjsToUse),x(subjsToUse));
            
            % plot reg line
            xlim = get(gca,'xlim');
            plot(xlim,[s.beta(2)*xlim(1)+s.beta(1) s.beta(2)*xlim(2)+s.beta(1)],'-k','linewidth',3)
            xlabel('Post-Stim Classifier Estimate % Change','fontsize',24)
            ylabel('Post-Stim Performance % Change','fontsize',24)
            grid on
            set(gca,'gridlinestyle',':')
            set(gca,'xlim',[xlim(1)-20 xlim(2)+20])
            
            % calc perm stats and plot perm reg line
            [r,p] = corr(x(subjsToUse),y(subjsToUse));
            r_perm = NaN(1,1000);
            p_perm = NaN(1,1000);
            permLine = NaN(1000,2);
            for i = 1:1000
                x_perm = xdata_perm{xDataType}(subjsToUse,i);%./postProbNonStim_perm(~bad,i);
                y_perm = ydata_perm{yDataType}(subjsToUse,i);%./nonStimPerf_perm(~bad,i);
                [r_perm(i),p_perm(i)] = corr(x_perm,y_perm);
                s = regstats(y_perm,x_perm,'linear',{'beta'});
                permLine(i,:) = s.beta;
            end
            
            goodBeta = sum(isnan(permLine) | isinf(permLine),2)==0;
            beta = median(permLine(goodBeta,:));            
            plot(xlim,[beta(2)*xlim(1)+beta(1) beta(2)*xlim(2)+beta(1)],'--k','linewidth',2)
            if filterInd < 3
            legend('EC','Hipp','MTL','Other','location','best')
            else
                legend('EC','Hipp','MTL','location','best')
            end
%             titleStr = sprintf('r(%d) = %.3f, p = %.3f, perm p = %.3f',sum(subjsToUse)-2,r,p,mean(p>p_perm));
            titleStr = sprintf('r(%d) = %.3f, p = %.3f, perm p = %.3f',sum(subjsToUse)-2,r,p,mean(r<r_perm));
            title(titleStr,'fontsize',18)
            set(gca,'TitleFontSizeMultiplier',1)
            set(gca,'TitleFontWeight','normal')
            set(gca,'fontsize',24)
            print(figs.(fname),'-depsc2','-loose')
        end
    end
end





texName = 'YC2_postStim_report.tex';
write_texfile(saveDir,texName,figs)

curr_dir = pwd;
cd(saveDir);
fprintf('Compiling pdf...\n');
unix(['pdflatex -shell-escape ' fullfile(saveDir, texName)]);
unix(['rm ' texName(1:end-3) 'aux']);
unix(['rm ' texName(1:end-3) 'log']);
fprintf('Done!\n');
cd(curr_dir);


% Start making the tex file
function write_texfile(saveDir,texName, figs)

% Write the document. If you do not have write permission, this will crash.
fid = fopen(fullfile(saveDir,texName),'w');

if fid==-1;
    error(sprintf('cannot open %s',texName))
end

% Write out the preamble to the tex doc. This is standard stuff and doesn't
% need to be changed
fprintf(fid,'\\documentclass[a4paper]{article} \n');
fprintf(fid,'\\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}\n');
fprintf(fid,'\\usepackage{graphicx,multirow} \n');
fprintf(fid,'\\usepackage{epstopdf} \n');
fprintf(fid,'\\usepackage[small,bf,it]{caption}\n');
fprintf(fid,'\\usepackage{subfig,amsmath} \n');
fprintf(fid,'\\usepackage{wrapfig} \n');
fprintf(fid,'\\usepackage{longtable} \n');
fprintf(fid,'\\usepackage{pdfpages}\n');
fprintf(fid,'\\usepackage{mathtools}\n');
fprintf(fid,'\\usepackage{array}\n');
fprintf(fid,'\\usepackage{enumitem}\n');
fprintf(fid,'\\usepackage{sidecap} \\usepackage{soul}\n');

% fprintf(fid,'\\setlength\\belowcaptionskip{5pt}\n');
fprintf(fid,'\n');
fprintf(fid,'\\addtolength{\\oddsidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\evensidemargin}{-.875in} \n');
fprintf(fid,'\\addtolength{\\textwidth}{1.75in} \n');
fprintf(fid,'\\addtolength{\\topmargin}{-.75in} \n');
fprintf(fid,'\\addtolength{\\textheight}{1.75in} \n');
fprintf(fid,'\n');
fprintf(fid,'\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}} \n');

fprintf(fid,'\\usepackage{fancyhdr}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
% fprintf(fid,'\\lhead{Report: %s }\n',strrep(subj,'_','\_'));
fprintf(fid,'\\rhead{YC2 Post Stimulation Report - Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
for s = 1:length(figs)
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_prob_perfScore_all);
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_prob_perfScore_good);
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_logOdds_perfScore_all);
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_logOdds_perfScore_good);
    fprintf(fid,'\\caption{Percent change in behavioral performance (\\textit{perfomance factor}) for stimulated items as a function of percent change in classifier output following stimulation. Each dot is a subject. Solid line is the true regression line, dashed line is the median permuation derived regression line. \\textbf{Top:} Classifier output measured as probability. \\textbf{Bottom:} Classifier output measured as log odds. \\textbf{Left:} All subjects. \\textbf{Right:} Subjects with significant YC1 decoding.}\n');
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_prob_perfScore_weights);        
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_logOdds_perfScore_weights);
    fprintf(fid,'\\caption{}\n');
    fprintf(fid,'\\end{figure}\n\n\n');    
    
    fprintf(fid,'\\clearpage\n\n\n');

    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_prob_RecChange_all);
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_prob_RecChange_good);
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_logOdds_RecChange_all);
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_logOdds_RecChange_good);
    fprintf(fid,'\\caption{Percent change in behavioral performance (\\textit{recall rate}) for stimulated items as a function of percent change in classifier output following stimulation. Each dot is a subject. Solid line is the true regression line, dashed line is the median permuation derived regression line. \\textbf{Top:} Classifier output measured as probability. \\textbf{Bottom:} Classifier output measured as log odds. \\textbf{Left:} All subjects. \\textbf{Right:} Subjects with significant YC1 decoding.}\n');
    fprintf(fid,'\\end{figure}\n\n\n');
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_prob_RecChange_weights);        
    fprintf(fid,'\\includegraphics[width=0.49\\textwidth]{%s}\n',figs(s).postStim_logOdds_RecChange_weights);
    fprintf(fid,'\\caption{}\n');
    fprintf(fid,'\\end{figure}\n\n\n');    
    
    fprintf(fid,'\\clearpage\n\n\n');    
    
end 

fprintf(fid,'\\end{document}\n\n\n');
fclose(fid);





