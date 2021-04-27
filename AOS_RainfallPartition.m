function [Runoff,Infl,Imax,NewCond] = AOS_RainfallPartition(P,Soil,FieldMngt,...
    InitCond,Slope_input,LUval,climateOS_Pori_day,extraP4cell_day)
% Function to partition rainfall into surface runoff and infiltration     
% using the curve number approach                                        

%% Store initial conditions for updating %%
NewCond = InitCond;

%% Calculate runoff %%
if (FieldMngt.Bunds == 0) || (FieldMngt.zBund < 0.001)
    % No bunds on field
    % Reset submerged days
    NewCond.DaySubmerged = 0;
    if Soil.AdjCN == 1 % Adjust CN for antecedent moisture
        % Check which compartment cover depth of top soil used to adjust
        % curve number
        comp_sto = find(Soil.Comp.dzsum>=Soil.zCN,1,'first');
        if isempty(comp_sto)
            comp_sto = Soil.nComp;
        end
        % Calculate weighting factors by compartment
        xx = 0;
        wrel = zeros(1,comp_sto);
        for ii = 1:comp_sto
            if Soil.Comp.dzsum(ii) > Soil.zCN
                Soil.Comp.dzsum(ii) = Soil.zCN;
            end
            wx = 1.016*(1-exp(-4.16*(Soil.Comp.dzsum(ii)/Soil.zCN)));
            wrel(ii) = wx-xx;
            if wrel(ii) < 0
                wrel(ii) = 0;
            elseif wrel(ii) > 1
                wrel(ii) = 1;
            end
            xx = wx;
        end
        % Calculate relative wetness of top soil
        wet_top = 0;
        for ii = 1:comp_sto
            layeri = Soil.Comp.Layer(ii);
            th = max(Soil.Layer.th_wp(layeri),InitCond.th(ii));
            wet_top = wet_top+(wrel(ii)*((th-Soil.Layer.th_wp(layeri))/...
                (Soil.Layer.th_fc(layeri)-Soil.Layer.th_wp(layeri))));
        end
        % Calculate adjusted curve number
        if wet_top > 1
            wet_top = 1;
        elseif wet_top < 0
            wet_top = 0;
        end
        CN = round(Soil.CNbot+(Soil.CNtop-Soil.CNbot)*wet_top);
    else % Curve number is not adjusted        
        CN = Soil.CN;
    end
    % Partition rainfall into runoff and infiltration (mm) 
    S = (25400/CN)-254;
    term = P-((5/100)*S);
    if term <= 0
        Runoff = 0;
        Infl = P;
    else
        Runoff = (term^2)/(P+(1-(5/100))*S);
        Infl = P-Runoff;
    end
else
    % Bunds on field, therefore no surface runoff 
    Runoff = 0;
    Infl = P;
end
    
%% Re-infiltration MVL 2017

if extraP4cell_day ~= 0
    % manning n
    minCC=0;
    maxCC=0.8;
    mina=0.03;
    maxb=0.3;
    if LUval == 1 % degraded LU
        if NewCond.CCadj ==0 % "bare soil"
            manning_n=0.03;
        else % crop development
            manning_n=(((maxb-mina)*(NewCond.CCadj-minCC))/(maxCC-minCC))+mina; % max CC = 0.8
        end
    else
        manning_n=0.4;
    end
    % manning V
    hydroR=0.005; % 5mm
    manningV=(Slope_input^(1/2)*hydroR^(2/3))/manning_n;
    
    % manning Time
    manningT=20/manningV;
    
    % Original Infiltration
    term_Pori = climateOS_Pori_day-((5/100)*S);
    if term_Pori <= 0
        Runoff_Pori = 0;
        Iori = P;
    else
        Runoff_Pori = (term_Pori^2)/(climateOS_Pori_day+(1-(5/100))*S);
        Iori = climateOS_Pori_day-Runoff_Pori;
    end
    
    % extra I
    Iextra=Infl-Iori;
    
    % Maximal Infiltration
    Imax=(Iextra/3600)*manningT;
    
    % new runoff
    Runoff=Runoff+(Infl-(Imax+Iori));
    
    % new Infiltration
    Infl=Imax+Iori;
else
    Imax=0;
end

end

