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

if any(strcmp(region,{'','all'}))
    fprintf('Using all electrodes for %s.\n',tal(1).subject)    
    return
end

% filter tal structure to just electrodes in a region of interest.
if ~isfield(tal,'locTag') || isempty([tal.locTag])
    fprintf('No loc tag information for %s.\n',tal(1).subject)
    if ~isempty(region)
        fprintf('Regional analysis requested by no localizations found, skipping %s.\n',tal(1).subject)
        tal = [];
        return
    end
else
    if strcmpi(region,'hipp')
        if any(~cellfun('isempty',regexpi({tal.locTag},['CA1|CA3|DG|sub'])))
            fprintf('Using only hippocampal electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['CA1|CA3|DG|sub'])));
        else
            fprintf('Using only hippocampal electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'ec')
        if any(~cellfun('isempty',regexpi({tal.locTag},['ec|erc'])))
            fprintf('Using only ec electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['ec|erc'])));
        else
            fprintf('Using only ec electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'CA1')
        if any(~cellfun('isempty',regexpi({tal.locTag},['CA1'])))
            fprintf('Using only CA1 electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['CA1'])));
        else
            fprintf('Using only CA1 electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'CA3')
        if any(~cellfun('isempty',regexpi({tal.locTag},['CA3'])))
            fprintf('Using only CA3 electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['CA3'])));
        else
            fprintf('Using only CA3 electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'DG')
        if any(~cellfun('isempty',regexpi({tal.locTag},['DG'])))
            fprintf('Using only DG electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['DG'])));
        else
            fprintf('Using only DG electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end               
    elseif strcmpi(region,'PHC')
        if any(~cellfun('isempty',regexpi({tal.locTag},['PHC'])))
            fprintf('Using only PHC electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['PHC'])));
        else
            fprintf('Using only PHC electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end                   
    elseif strcmpi(region,'PRC')
        if any(~cellfun('isempty',regexpi({tal.locTag},['PRC'])))
            fprintf('Using only PRC electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['PRC'])));
        else
            fprintf('Using only PRC electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end      
    elseif strcmpi(region,'Sub')
        if any(~cellfun('isempty',regexpi({tal.locTag},['Sub'])))
            fprintf('Using only Sub electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['Sub'])));
        else
            fprintf('Using only Sub electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'mtl')
        if any(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
            fprintf('Using only mtl electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])));
        else
            fprintf('Using only mtl electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    else
        if any(~cellfun('isempty',regexpi({tal.locTag},[region])))
            fprintf('Using only %s electrodes for %s.\n',region,tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi({tal.locTag},[region])));
        else
            fprintf('Using only %s electrodes for %s...NONE FOUND.\n',region,tal(1).subject)
            tal = [];
            return
        end
    end
end