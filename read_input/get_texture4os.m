function texture=get_texture4os(lithology,jj,ii,xx1,yy1)
%sand - clay - OM
% col=[0.4,0.2,3]; % 1=colluvial
% lime=[0.3,0.2,2]; % 2=Limestone
% oph=[0.7,0.09,2]; % 3=Ophiolite
% con=[0.5,0.15,2]; % 4=Conglomerate
col=[0.4,0.2,3]; % 1=colluvial
lime=[0.3,0.2,2]; % 2=Limestone
con_oph=[(0.7+0.5)/2,(0.09+0.15)/2,2]; % 3=Ophiolite
marl=[0.15,0.28,2]; % 4=Conglomerate
textureMAT=[col;lime;con_oph;marl];
%select sand - clay - OM based on lithology
texture=textureMAT(lithology(jj+yy1-1,ii+xx1-1),:);

%SCL - comparison Batch
% texture=[0.08,0.30,2];
end

