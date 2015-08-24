function tal = filterTalByRegion(tal,region)
% function tal = filterTalByRegion(region)
%
% Inputs:   region - string of region of interest
%                    'hipp' includes CA1|CA3|DG|sub
%                    'ec' includes ec|erc 
%                    'mtl' includes HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc
%
% Output: filtered tal structure
%
% TO DO: add support for more regions

% filter tal structure to just electrodes in a region of interest.
if ~isfield(tal,'locTag')
    fprintf('No loc tag information for %s.\n',subj)
    if ~isempty(region)
        fprintf('Regional analysis requested by no localizations found, skipping %s.\n',subj)
        return
    end
else
    if strcmpi(region,'hipp')
        if any(~cellfun('isempty',regexpi({tal.locTag},['CA1|CA3|DG|sub'])))
            fprintf('Using only hippocampal electrodes for %s.\n',subj)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['CA1|CA3|DG|sub'])));
        else
            fprintf('Using only hippocampal electrodes for %s...NONE FOUND.\n',subj)
            return
        end
    elseif strcmpi(region,'ec')
        if any(~cellfun('isempty',regexpi({tal.locTag},['ec|erc'])))
            fprintf('Using only ec electrodes for %s.\n',subj)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['ec|erc'])));
        else
            fprintf('Using only ec electrodes for %s...NONE FOUND.\n',subj)
            return
        end
    elseif strcmpi(region,'mtl')
        if any(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
            fprintf('Using only mtl electrodes for %s.\n',subj)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])));
        else
            fprintf('Using only mtl electrodes for %s...NONE FOUND.\n',subj)
            return
        end
    else
       fprintf('Region string not recognized.\n')
    end    
end