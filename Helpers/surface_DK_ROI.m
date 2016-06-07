function anat = surface_DK_ROI(subjTalEvents,subj)
% updated 12.12.13

c_temp = {subjTalEvents.channel};
chans = vertcat(c_temp{:});

IFG = zeros(1,length(chans));
MFG = zeros(1,length(chans));
SFG = zeros(1,length(chans));

MTL = zeros(1,length(chans));
Hipp = zeros(1,length(chans));
TC = zeros(1,length(chans));

IPC = zeros(1,length(chans));
SPC = zeros(1,length(chans));

OC = zeros(1,length(chans));

subjTalEv = [subjTalEvents(:).avgSurf];
%subjElec   = [subjTalEv.channel]';
subjXYZ    = [[subjTalEv.x]' [subjTalEv.y]' [subjTalEv.z]'];

subjDK = {subjTalEv.anatRegion}';

subjTalEv = [subjTalEvents(:)];

if isfield(subjTalEv,'locTag')
    subjlocTag = {subjTalEv.locTag};
    
    Hipp_locTag = {'Left CA1','Left CA2','Left CA3','Left DG','Left Sub',...
        'Right CA1','Right CA2','Right CA3','Right DG','Right Sub'};
    
    MTL_locTag = {'Left PRC','Right PRC','Right EC','Right PHC','Left EC','Left PHC'};
else
    hipp_elec = get_mtl_leads(subj,'hipp',0);
    phc_elec = get_mtl_leads(subj,'phc',0);
end


[a,b] = size(subjXYZ);
if a > 0 && b > 0
    LEFT = logical(subjXYZ(:,1) < 0);
    RIGHT = logical(subjXYZ(:,1) > 0);
else
    LEFT = [];
    RIGHT = [];
end


%%%%%%%%%%%%%%%%%% BOTH BIPOLAR ELECTRODES HAVE TO BE IN THE ROI.
% not for hipp but for cortical
for c = 1:length(chans)
    dk = subjDK(c);
    
    if isfield(subjTalEv,'locTag')
        locTag = subjlocTag(c);
        if isempty(locTag{1})
            locTag = {''};
        end
    else
        first = chans(c,1);
        second = chans(c,2);
    end
    
    if ~isempty(dk)
              
%         aPFC(c) = strcmp('Brodmann area 10',dk) | strcmp('Brodmann area 11',dk);
        IFG(c) = strcmp('parsopercularis',dk) | strcmp('parsorbitalis',dk) |...
            strcmp('parstriangularis',dk);
        
        MFG(c) = strcmp('caudalmiddlefrontal',dk) | strcmp('rostralmiddlefrontal',dk);
        
        SFG(c) = strcmp('superiorfrontal',dk);
        
        TC(c) = strcmp('middletemporal',dk) | strcmp('inferiortemporal',dk);% |...
            %strcmp('Brodmann area 37',dk);
        
        IPC(c) = strcmp('inferiorparietal',dk) | strcmp('supramarginal',dk);
        
        SPC(c) = strcmp('superiorparietal',dk) | strcmp('precuneus',dk);
        
        OC(c) = strcmp('lateraloccipital',dk) | strcmp('lingual',dk) |...
            strcmp('cuneus',dk) | strcmp('pericalcarine',dk);
        
        if isfield(subjTalEv,'locTag')
            MTL(c) = ismember(locTag,MTL_locTag);
            Hipp(c) = ismember(locTag,Hipp_locTag);
        else
            % Hipp
            for i = 1:length(hipp_elec)
                for j = 1:length(hipp_elec)
                    if first == hipp_elec(i)|| second == hipp_elec(j)
                        Hipp(c) = 1;
                    end
                end
            end
            
            % MTL
            if isempty(phc_elec)
                MTL(c) =  ((strcmp('parahippocampal',dk)) || (strcmp('entorhinal',dk)))  & ~Hipp(c);                
            else
                for i = 1:length(phc_elec)
                    for j = 1:length(phc_elec)
                        if first == phc_elec(i)|| second == phc_elec(j)
                            MTL(c) = 1;
                        end
                    end                    
                end
            end
            
            %MTL(c) =  PHC(c) || Hipp(c);
            MTL(c) = ((strcmp('parahippocampal',dk)) || (strcmp('entorhinal',dk))) & ~Hipp(c);

        end
        
    end
end

IFG(MTL|Hipp)=0;
MFG(MTL|Hipp)=0;
SFG(MTL|Hipp)=0;
TC(MTL|Hipp)=0;
IPC(MTL|Hipp)=0;
SPC(MTL|Hipp)=0;
OC(MTL|Hipp)=0;

anat.ROI.IFG = IFG;
anat.ROI.MFG = MFG;
anat.ROI.SFG = SFG;
anat.ROI.MTL = MTL;
anat.ROI.Hipp = Hipp;
anat.ROI.TC = TC;
anat.ROI.IPC = IPC;
anat.ROI.SPC = SPC;
anat.ROI.OC = OC;

anat.ROI.remain = ~(IFG | MFG | SFG | MTL | Hipp | TC | IPC | SPC | OC);

anat.hemi.right = RIGHT';
anat.hemi.left = LEFT';
