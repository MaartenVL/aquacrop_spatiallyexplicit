function [ParamStruct] = readCrop4os(AQpath)
load([AQpath 'ParamOS.mat']); % ParamStruct
%change GDD info, planting and harvest dates 
PlantingDate='22/09';
HarvestDate='07/07';
% Emergence=435;
% MaxCC=6293;
% MaxRD=5713;
% StartCanSen=7192;
% Maturity=7830;
% Flowering=6670;
% BuildupHI=348;
% DurFlow=348;
% Emergence=289;
% MaxCC=800;
% MaxRD=750;
% StartCanSen=2835;
% Maturity=3390;
% DurFlowering=264;
% time2YieldForm=1073;
% DurYieldForm=300;

% fprintf( fid, '%s\n','   100         : Calendar Days: from sowing to emergence');
% fprintf( fid, '%s\n','   217         : Calendar Days: from sowing to maximum rooting depth');
% fprintf( fid, '%s\n','   222         : Calendar Days: from sowing to start senescence');
% fprintf( fid, '%s\n','   300         : Calendar Days: from sowing to maturity (length of crop cycle)');
% fprintf( fid, '%s\n','   200         : Calendar Days: from sowing to flowering');
% fprintf( fid, '%s\n','    90         : Length of the flowering stage (days)');
% fprintf( fid, '%s\n','     1         : Crop determinancy linked with flowering');
% fprintf( fid, '%s\n','   100         : Excess of potential fruits (%)');
% fprintf( fid, '%s\n','    27         : Building up of Harvest Index starting at flowering (days)');

%From field work
% emergence=14;
% flowering=233;
% Maturity=276;

Emergence=14;
MaxCC=204;
MaxRD=240;
StartCanSen=250;
Maturity=276;
DurFlowering=20;
time2YieldForm=233; % idem time to flowering
DurYieldForm=37;


%fill in Crop File - winter Barley from AQbatch5 - some info is missing in
%the AQ v5.0 generated CropFile, I leave this as provided in the AQ OS
%cropfile, might be wrong for a barley crop...
ParamStruct.Crop.Maize.CalendarType=1;
ParamStruct.Crop.Maize.SwitchGDD=1;
ParamStruct.Crop.Maize.PlantingDate=PlantingDate;
ParamStruct.Crop.Maize.HarvestDate=HarvestDate;
ParamStruct.Crop.Maize.Emergence=Emergence;
ParamStruct.Crop.Maize.MaxRooting=MaxRD;
ParamStruct.Crop.Maize.Senescence=StartCanSen;
ParamStruct.Crop.Maize.Maturity=Maturity;
ParamStruct.Crop.Maize.HIstart=time2YieldForm;
ParamStruct.Crop.Maize.Flowering=DurFlowering;
ParamStruct.Crop.Maize.YldForm=DurYieldForm;

ParamStruct.Crop.Maize.fshape_w3=3;
ParamStruct.Crop.Maize.fshape_w2=3; 
ParamStruct.Crop.Maize.fshape_w1=3;
ParamStruct.Crop.Maize.p_lo1=0.65;
ParamStruct.Crop.Maize.p_up1=0.2;
ParamStruct.Crop.Maize.p_up2=0.6;
ParamStruct.Crop.Maize.p_up3=0.55;
ParamStruct.Crop.Maize.p_up4 =0.85;
ParamStruct.Crop.Maize.exc=100;
ParamStruct.Crop.Maize.dHI_pre=5;
ParamStruct.Crop.Maize.HI0=0.33;
ParamStruct.Crop.Maize.b_HI=5;
ParamStruct.Crop.Maize.a_HI=10;
ParamStruct.Crop.Maize.WP=15;
ParamStruct.Crop.Maize.fage=0.150;
ParamStruct.Crop.Maize.Kcb=1.1;
ParamStruct.Crop.Maize.CGC=0.03209; 
ParamStruct.Crop.Maize.CDC=0.07697;
ParamStruct.Crop.Maize.CCx=0.8;
ParamStruct.Crop.Maize.PlantPop=1500000;
ParamStruct.Crop.Maize.SeedSize=1.5;
ParamStruct.Crop.Maize.SxBotQ=0.006;
ParamStruct.Crop.Maize.SxTopQ=0.019; 
ParamStruct.Crop.Maize.Zmax=1.3;
ParamStruct.Crop.Maize.GDD_up=14;
ParamStruct.Crop.Maize.Tmin_lo=0; %??? no idea
ParamStruct.Crop.Maize.Tmin_up=5;
ParamStruct.Crop.Maize.Tmax_up=35;
ParamStruct.Crop.Maize.Tbase=0;
ParamStruct.Crop.Maize.Tupp=15;
ParamStruct.Crop.Maize.Aer=15;

%aanmaken van cropChoices - werk dit?
ParamStruct.Crop.Maize.CropChoices={'Maize'};

ParamStruct;

end

