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
    locTags          = {tal.locTag};
    missing          = cellfun('isempty',locTags);
    locTags          = strcat(locTags,'_mtl');
    lobes            = {tal.Loc2};    
    badLobes         = cellfun('isempty',regexpi(lobes,['temporal|occipital|parietal|frontal|limbic']));    
    lobes(badLobes)  = repmat({''},1,sum(badLobes));
    locTags(missing) = lobes(missing);
    
    if strcmpi(region,'hipp')
        if any(~cellfun('isempty',regexpi(locTags,['CA1|CA3|DG|sub'])))
            fprintf('Using only hippocampal electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['CA1|CA3|DG|sub'])));
        else
            fprintf('Using only hippocampal electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'ec')
        if any(~cellfun('isempty',regexpi(locTags,['ec|erc'])))
            fprintf('Using only ec electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['ec|erc'])));
        else
            fprintf('Using only ec electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'CA1')
        if any(~cellfun('isempty',regexpi(locTags,['CA1'])))
            fprintf('Using only CA1 electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['CA1'])));
        else
            fprintf('Using only CA1 electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'CA3')
        if any(~cellfun('isempty',regexpi(locTags,['CA3'])))
            fprintf('Using only CA3 electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['CA3'])));
        else
            fprintf('Using only CA3 electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'DG')
        if any(~cellfun('isempty',regexpi(locTags,['DG'])))
            fprintf('Using only DG electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['DG'])));
        else
            fprintf('Using only DG electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end               
    elseif strcmpi(region,'PHC')
        if any(~cellfun('isempty',regexpi(locTags,['PHC'])))
            fprintf('Using only PHC electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['PHC'])));
        else
            fprintf('Using only PHC electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end                   
    elseif strcmpi(region,'PRC')
        if any(~cellfun('isempty',regexpi(locTags,['PRC'])))
            fprintf('Using only PRC electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['PRC'])));
        else
            fprintf('Using only PRC electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end      
    elseif strcmpi(region,'Sub')
        if any(~cellfun('isempty',regexpi(locTags,['Sub'])))
            fprintf('Using only Sub electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['Sub'])));
        else
            fprintf('Using only Sub electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'mtl')
        if any(~cellfun('isempty',regexpi(locTags,['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
            fprintf('Using only mtl electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])));
        else
            fprintf('Using only mtl electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'frontal')        
        if any(~cellfun('isempty',regexpi(locTags,['Frontal'])))
            fprintf('Using only frontal electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['Frontal'])));
        else
            fprintf('Using only frontal electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end
    elseif strcmpi(region,'parietal')        
        if any(~cellfun('isempty',regexpi(locTags,['parietal'])))
            fprintf('Using only parietal electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['parietal'])));
        else
            fprintf('Using only parietal electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end        
    elseif strcmpi(region,'occipital')        
        if any(~cellfun('isempty',regexpi(locTags,['Occipital'])))
            fprintf('Using only Occipital electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['Occipital'])));
        else
            fprintf('Using only Occipital electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end            
    elseif strcmpi(region,'temporal')        
        if any(~cellfun('isempty',regexpi(locTags,['Temporal'])))
            fprintf('Using only Temporal electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['Temporal'])));
        else
            fprintf('Using only Temporal electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end         
%     elseif strcmpi(region,'temporal_nonmtl')        
%         if any(~cellfun('isempty',regexpi(locTags,['Temporal'])) & cellfun('isempty',regexpi(locTags,['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
%             fprintf('Using only Temporal non-MTL electrodes for %s.\n',tal(1).subject)
%             tal = tal(~cellfun('isempty',regexpi(locTags,['Temporal'])) & cellfun('isempty',regexpi(locTags,['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])));
%         else
%             fprintf('Using only Temporal Non-MTL electrodes for %s...NONE FOUND.\n',tal(1).subject)
%             tal = [];
%             return
%         end         
    elseif strcmpi(region,'limbic')        
        if any(~cellfun('isempty',regexpi(locTags,['limbic'])))
            fprintf('Using only limbic electrodes for %s.\n',tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,['limbic'])));
        else
            fprintf('Using only limbic electrodes for %s...NONE FOUND.\n',tal(1).subject)
            tal = [];
            return
        end       
%     elseif strcmpi(region,'limbic_nonmtl')        
%         if any(~cellfun('isempty',regexpi(locTags,['limbic'])) & cellfun('isempty',regexpi(locTags,['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])))
%             fprintf('Using only limbic non-MTL electrodes for %s.\n',tal(1).subject)
%             tal = tal(~cellfun('isempty',regexpi(locTags,['limbic'])) & cellfun('isempty',regexpi(locTags,['HC|ec|hipp|CA1|CA3|DG|sub|amy|phc|prc|BA36|erc'])));
%         else
%             fprintf('Using only limbic Non-MTL electrodes for %s...NONE FOUND.\n',tal(1).subject)
%             tal = [];
%             return
%         end          
    else        
        if any(~cellfun('isempty',regexpi(locTags,[region])))
            fprintf('Using only %s electrodes for %s.\n',region,tal(1).subject)
            tal = tal(~cellfun('isempty',regexpi(locTags,[region])));
        else
            fprintf('Using only %s electrodes for %s...NONE FOUND.\n',region,tal(1).subject)
            tal = [];
            return
        end
    end
end