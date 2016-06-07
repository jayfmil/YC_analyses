function res = YC1_DifficultyPerformanceReg(subjs,task,saveDir)

% get list of YC subjects if non given
if ~exist('subjs','var') || isempty(subjs)
    subjs = get_subs(task);
end
subjs_all = get_subs('RAM_YC1');

saveDir = fullfile(saveDir,task);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

figDir = fullfile(saveDir,'figs');
if ~exist(figDir,'dir')
    mkdir(figDir);
end

figs = [];
for s = 1:length(subjs)
    
    try
        % load events for subject
        events = get_sub_events(task,subjs{s});
        events = addExtraYCFields(events);
        
        % behavioral performance (1 = best, 0 = worst)
        pf = 1-[events(strcmp({events.type},'NAV_TEST')).testError];
        
        testEvents = find(strcmp({events.type},'NAV_TEST'));
        if length(testEvents) < 24
            continue
        end
        
        %%%% test metrics that may contribute to trial difficulty %%%
        % learnOneTestDist  - distance between learning trial one and test start
        % learnTwoTestDist  - distance between learning trial two and test start
        % learnOneTwoDist   - distance between learning trial one and two start
        % testObjectDist    - distance between test start and object location
        % learnOneTestAngle - cosine of angle between learning one drive and test
        % learnTwoTestAngle - cosine of angle between learning two drive and test
        
        % add in angle between test start and correct location
        % also some measure of the object location itself? I thought about
        % using the pf of all other subjects for that location, or
        % something like the rank of that location relative to other
        % locations, but that doesn't feel right... what about a measure of
        % distance of object to center? though maybe it is ok the first
        % way..
        
        
        % mean distance to two closest walls
        % do for each session sep
        
        % load errors for all other subjs
        other_subjs = subjs_all(~strcmp(subjs_all,subjs{s}));
        %[allErrors,allObjectLocs,allEucErrors,subjVec] = YC1_loadAllSubjErrors(0,other_subjs);
        
        
        learnOneTestDist  = NaN(1,length(testEvents));
        learnTwoTestDist  = NaN(1,length(testEvents));
        learnOneTwoDist   = NaN(1,length(testEvents));
        testObjectDist    = NaN(1,length(testEvents));
        learnOneTestAngle = NaN(1,length(testEvents));
        learnTwoTestAngle = NaN(1,length(testEvents));
        testObjectAngle   = NaN(1,length(testEvents));
        %objLocDiff        = NaN(1,length(testEvents));
        objWallDist       = NaN(1,length(testEvents));
        sessions          = NaN(1,length(testEvents));
        for i = 1:length(testEvents)
            learnLoc1     = events(testEvents(i)-2).startLocs;
            learnLoc2     = events(testEvents(i)-1).startLocs;
            testLoc       = events(testEvents(i)).startLocs;
            objLoc        = events(testEvents(i)).objLocs;
            sessions(i)   = events(testEvents(i)).session;
            
            % euc distance between learning trial one and and test
            learnOneTestDist(i) = sqrt(sum((learnLoc1 - testLoc).^2));
            %         learnOneTestDist(i) = calc_YC_error(learnLoc1,testLoc);
            
            % euc distance between learning trial two and and test
            learnTwoTestDist(i) = sqrt(sum((learnLoc2 - testLoc).^2));
            %         learnTwoTestDist(i) = calc_YC_error(learnLoc2,testLoc);
            
            % euc distance between learning trial one and two starts
            learnOneTwoDist(i) = sqrt(sum((learnLoc2 - learnLoc1).^2));
            %         learnOneTwoDist(i) = calc_YC_error(learnLoc2,learnLoc1);
            
            % euc distance between test start and correct location
            testObjectDist(i) = sqrt(sum((objLoc - testLoc).^2));
            %         testObjectDist(i) = calc_YC_error(objLoc,testLoc);
            
            % angle measures (sometimes Path is not filled in?)
            if ~isempty(events(testEvents(i)).Path.direction)
                learnOneAngle = mod(events(testEvents(i)-2).Path.direction(end),360);
                learnTwoAngle = mod(events(testEvents(i)-1).Path.direction(end),360);
                testAngle     = mod(events(testEvents(i)).Path.direction(1),360);
                learnOneTestAngle(i) = -1*cosd(learnOneAngle - testAngle);
                learnTwoTestAngle(i) = -1*cosd(learnTwoAngle - testAngle);
                
                % figure out anlge to correct location. 
                deltaLoc = objLoc - testLoc;
                ang = (events(testEvents(i)).Path.direction(1) - mod(atan2d(deltaLoc(2),deltaLoc(1))+180,180));
                ang = abs(90-mod(ang + 180,180));
                testObjectAngle(i) = ang;                
            end
            
            % object location
            %x_pos = objLoc(1);
            %y_pos = objLoc(2);
            %near = sqrt((allObjectLocs(:,1) - x_pos).^2 + (allObjectLocs(:,2) - y_pos).^2) < 5;
            %objLocDiff(i) = mean(allErrors(near));            
            
            % object locaiton distance to two closest walls
            wallDists = min(abs([objLoc(1) - 32.4 objLoc(1) + 32.4; objLoc(2) - 18 objLoc(2) + 18]'));
            objWallDist(i) = mean(wallDists);
            
        end
        
%         % create predictors matrix and normalize between zero and one
%         x = [learnOneTestDist' learnTwoTestDist' learnOneTwoDist' testObjectDist' learnOneTestAngle' learnTwoTestAngle' objWallDist' testObjectAngle'];
%         x = (x - repmat(min(x),size(x,1),1))./repmat(max(x - repmat(min(x),size(x,1),1)),size(x,1),1);
%         
        % response variable
        y = pf';
        
        % use matlab linear model class to fit model
        %tbl = table(x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6),x(:,7),x(:,8),y,'VariableNames',{'Learn1Dist','Learn2Dist','LearnDist','testDist','Learn1Angle','Learn2Angle','ObjectWallDist','testObjectAngle','PF'});
        %lm = fitlm(tbl,'PF~Learn1Dist+Learn2Dist+LearnDist+testDist+Learn1Angle+Learn2Angle+testObjectAngle+ObjectWallDist+testDist:testObjectAngle+Learn1Dist:Learn1Angle+Learn2Dist:Learn2Angle')
        
        % now predict the responses using hold out trials, looping over each
        % trial
        predictedPF = NaN(1,length(testEvents));
        resids     = NaN(1,length(testEvents));
        
        % add new field to events
        [events.predictedPF] = deal(NaN);
        [events.residual] = deal(NaN);

        
        uniqSessions = unique(sessions);
        rs = NaN(1,length(unique(sessions)));
        for sess = 1:length(unique(sessions))
            testIndsSess   = sessions == uniqSessions(sess);
            testEventsSess = testEvents(testIndsSess);
            
            
            % create predictors matrix and normalize between zero and one
            x = [learnOneTestDist(testIndsSess)' learnTwoTestDist(testIndsSess)' ...
                 learnOneTwoDist(testIndsSess)' testObjectDist(testIndsSess)' ...
                 learnOneTestAngle(testIndsSess)' learnTwoTestAngle(testIndsSess)' ...
                 objWallDist(testIndsSess)' testObjectAngle(testIndsSess)'];
            x = (x - repmat(min(x),size(x,1),1))./repmat(max(x - repmat(min(x),size(x,1),1)),size(x,1),1);
            
            % response variable
            ySess = pf(testIndsSess)';
            
            % use matlab linear model class to fit model
            tbl = table(x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6),x(:,7),x(:,8),ySess,'VariableNames',{'Learn1Dist','Learn2Dist','LearnDist','testDist','Learn1Angle','Learn2Angle','ObjectWallDist','testObjectAngle','PF'});
            lm = fitlm(tbl,'PF~Learn1Dist+Learn2Dist+LearnDist+testDist+Learn1Angle+Learn2Angle+testObjectAngle+ObjectWallDist+testDist:testObjectAngle+Learn1Dist:Learn1Angle+Learn2Dist:Learn2Angle')
            rs(sess) = lm.Rsquared.Ordinary;
            
            
            
            for trial = 1:length(testEventsSess)
                
                % training data set
                % trainTrials = setdiff(1:length(testEventsSess),trial);
                trainTrials   = find(testIndsSess);
                testTrial     = trainTrials(trial);
                trainTrials(trial) = [];
                
                xTrain = [learnOneTestDist(trainTrials)' learnTwoTestDist(trainTrials)' ...
                    learnOneTwoDist(trainTrials)' testObjectDist(trainTrials)' ...
                    learnOneTestAngle(trainTrials)' learnTwoTestAngle(trainTrials)' ...
                    objWallDist(trainTrials)' testObjectAngle(trainTrials)'];
                
                xMin   = min(xTrain);
                xMax   = max(xTrain - repmat(xMin,size(xTrain,1),1));
                xTrain = bsxfun(@rdivide,bsxfun(@minus,xTrain,xMin),xMax);
                
                % create model
                tbl = table(xTrain(:,1),xTrain(:,2),xTrain(:,3),xTrain(:,4),xTrain(:,5),xTrain(:,6),xTrain(:,7),xTrain(:,8),y(trainTrials),'VariableNames',{'Learn1Dist','Learn2Dist','LearnDist','testDist','Learn1Angle','Learn2Angle','ObjectWallDist','testObjectAngle','PF'});
                lmTrain = fitlm(tbl,'PF~Learn1Dist+Learn2Dist+LearnDist+testDist+Learn1Angle+Learn2Angle+testObjectAngle+ObjectWallDist+testDist:testObjectAngle+Learn1Dist:Learn1Angle+Learn2Dist:Learn2Angle');
                
                % predict
                xTest = [learnOneTestDist(testTrial)' learnTwoTestDist(testTrial)' learnOneTwoDist(testTrial)' testObjectDist(testTrial)' learnOneTestAngle(testTrial)' learnTwoTestAngle(testTrial)' objWallDist(testTrial)' testObjectAngle(testTrial')];
                xTest = (xTest - xMin)./xMax;
                predictedPF(testTrial) = predict(lmTrain,xTest);
                
                [events(testEventsSess(trial)-2:testEventsSess(trial)).predictedPF] = deal(predictedPF(testTrial));
                [events(testEventsSess(trial)-2:testEventsSess(trial)).residual] = deal(y(trial) -  predictedPF(testTrial));
                
            end
        end
        % save out events with added fields
        ev_file = fullfile(saveDir,[subjs{s} '_events.mat']);
        save(ev_file,'events');
        
        
        % make and save figures
        figs_subj = [];
        figs_subj.subj = subjs{s};
        
        
        figure(1)
        clf
        hold on
        titleStr = [strrep(subjs{s},'_',' ') ': '];
        for sess = 1:length(uniqSessions)
            if sess < length(uniqSessions)
                titleStr = sprintf('%s R2-%d=%.3f, ',titleStr,sess,rs(sess));
            else
                titleStr = sprintf('%s R2-%d=%.3f',titleStr,sess,rs(sess));
            end
            scatter(y(sessions==uniqSessions(sess)),predictedPF(sessions==uniqSessions(sess))',80,'filled')
        end
        set(gca,'fontsize',20)
        set(gca,'xlim',[0 1]);
        set(gca,'ylim',[0 1]);
        xlabel('Performance Factor','fontsize',20)
        ylabel('Predicted Performance Factor','fontsize',20)
        grid on
        set(gca,'gridlinestyle',':')
        %titleStr = sprintf('%s: rsquared = %.3f, adj. rsquared = %.3f',strrep(subjs{s},'_',' '),lm.Rsquared.Ordinary,lm.Rsquared.Adjusted);
        h2=title(titleStr);
        h2.FontWeight = 'normal';
        h2.FontSize = 17;
        legend( cellfun(@num2str, num2cell(1:length(uniqSessions)), 'UniformOutput', false),'location','eastoutside')
        fname = fullfile(figDir,[subjs{s} '_pf_x_pfHat.eps']);
        figs_subj.pfHat = fname;
        print('-depsc2','-loose',figs_subj.pfHat)
        
        figure(2)
        clf
        resids = y - predictedPF';
        %scatter(resids,y,80,'filled')
        matched = (resids < 0 & y < median(y)) | (resids > 0 & y > median(y));
        scatter(resids(matched),y(matched),80,'filled','k')
        hold on
        h=scatter(resids(~matched),y(~matched),80,'filled','r');
        h.CData = [.8 .2 .2];
        ylim = get(gca,'ylim');
        xlim = get(gca,'xlim');
        plot(xlim,[median(y) median(y)],'--k','linewidth',2)
        plot([0 0],ylim,'--k','linewidth',2)
        set(gca,'fontsize',20)
        set(gca,'xlim',xlim);
        set(gca,'ylim',ylim);
        xlabel('Residual','fontsize',20)
        ylabel('Performance Factor','fontsize',20)
        grid on
        set(gca,'gridlinestyle',':')
        titleStr = sprintf('Percent overlap: %.3f',mean(matched)*100);
        h2=title(titleStr);
        h2.FontWeight = 'Normal';
        fname = fullfile(figDir,[subjs{s} '_resid.eps']);
        figs_subj.resid = fname;
        print('-depsc2','-loose',figs_subj.resid)
        
        figure(3)
        clf
        subplot(2,1,1)
        [n,x] = hist(y,20);
        bar(x,n./sum(n),1,'w','linewidth',2)
        h=title('PF Dist.');
        h.FontWeight = 'Normal';
        grid on
        set(gca,'gridlinestyle',':')
        ylabel('Prob.','fontsize',20)
        set(gca,'fontsize',20)
        subplot(2,1,2)
        [n,x] = hist(resids,20);
        bar(x,n./sum(n),1,'w','linewidth',2)
        h=title('Resid Dist.');
        h.FontWeight = 'Normal';
        ylabel('Prob.','fontsize',20)
        grid on
        set(gca,'gridlinestyle',':')
        set(gca,'fontsize',20)
        fname = fullfile(figDir,[subjs{s} '_hist.eps']);
        figs_subj.hist = fname;
        print('-depsc2','-loose',figs_subj.hist)  
        
        figs = [figs figs_subj];
        
    end
end


texName = [task '_residReport.tex'];
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
fprintf(fid,'\\rhead{Date created: %s}\n',date);

fprintf(fid,'\\usepackage{hyperref}\n');

% Start the document
fprintf(fid,'\\begin{document}\n\n\n');

% fprintf(fid,'\\hypertarget{%s}{}\n',region{1});

% This section writes the figures
for s = 1:length(figs)
    
    fprintf(fid,'\\begin{figure}[!h]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\includegraphics[width=0.32\\textwidth]{%s}\n',figs(s).pfHat);
    fprintf(fid,'\\includegraphics[width=0.32\\textwidth]{%s}\n',figs(s).resid);
    fprintf(fid,'\\includegraphics[width=0.32\\textwidth]{%s}\n',figs(s).hist);
    fprintf(fid,'\\caption{%s}\n\n',strrep(figs(s).subj,'_',' '));
    fprintf(fid,'\\end{figure}\n\n\n');
    if mod(s,4) == 0
        fprintf(fid,'\\clearpage\n\n\n');
    end
end

fprintf(fid,'\\end{document}\n\n\n');









