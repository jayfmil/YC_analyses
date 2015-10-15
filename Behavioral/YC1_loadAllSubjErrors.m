function [allErrors,allObjectLocs,allEucErrors,subjVec] = YC1_loadAllSubjErrors(recalcError)
% function [allErrors,allObjectLocs,allEucErrors,subjVec] = YC1_loadAllSubjErrors(recalcError)
%
% This functions loops over all the YC1 subjects in the dataset and returns
% a vector of all the errors made and a n x 2 matrix of all the
% corresponding object locations.

if ~exist('recalcError','var') || isempty(recalcError)
    recalcError = 0;
end

% get list of YC subjects
subjs = get_subs('RAM_YC1');

% loop over each subject
errors = cell(1,length(subjs));
objLocs = cell(1,length(subjs));
eucErrors = cell(1,length(subjs));

disp('Loading all errors.')
for s = 1:length(subjs);    
    % process subject    
    [errors{s},eucErrors{s},objLocs{s}] = process_subj(subjs{s},recalcError);    
end

allErrors     = vertcat(errors{:});
allEucErrors  = vertcat(eucErrors{:});
allObjectLocs = vertcat(objLocs{:});

% make a vector to keep track of which entries are from which subject
% number
countsPerSubj = cellfun(@numel,errors);
subjVec       = zeros(1,sum(countsPerSubj));
subjVec(cumsum(countsPerSubj) - countsPerSubj + 1) = 1;
subjVec = cumsum(subjVec);


if length(allErrors) ~= length(allObjectLocs)
    keyboard
end


%errs = NaN(61,31);
%for x = -30:30
%    for y = -15:15        
%        near = sqrt((allObjectLocs(:,1) - x).^2 + (allObjectLocs(:,2) ...
%                                                   - y).^2) < 10;
%        errs(x+31,y+16) = mean(allErrors(near));        
%    end
%end
%keyboard


function [errors,eucErrors,objLocs] = process_subj(subj,doRecalc)


% load events and filter to test trials
events = get_sub_events('RAM_YC1',subj);
events = events(strcmp({events.type},'NAV_TEST'));

% get all object locations and response error
objLocs = vertcat(events.objLocs);
respLocs = vertcat(events.respLocs);
eucErrors = [events.respDistErr]';

% calculate burgess style error metric based on distribution of possible
% errors
if doRecalc
    errors = NaN(length(events),1);
    for e = 1:length(events)
        errors(e) = calc_YC_error(objLocs(e,:),respLocs(e,:));
    end
else    
    errors = [events.respPerformanceFactor]';
end




