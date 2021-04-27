function  write_mvl_input( climateOS,dates_mvl,SoilInfo,CO2atTime,SWCtemp,AQpath,j,w,CropReadFile,CropReadFileTrees,LUval)
path=[AQpath 'Input/'];


%% Crop File
if LUval == 1 % degraded LU
    [fileCrop]=create_crop(CropReadFile,path,j,w);
else %undegraded LU
    CropReadFile.CropRead=CropReadFileTrees.CropwheatCal4Paleo4trees;
    [fileCrop]=create_crop(CropReadFile,path,j,w);
end

%% Crop Rotation
%wordt niet gelezen dus hier ook niet aanmaken
%% Climate
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['Weather_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileClim=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Weather input time-series for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Day	Month	Year	MinTemp	MaxTemp	Precipitation	ReferenceET %%');
for i=1:365
    fprintf( fid, '%i\t%i\t%i\t%f\t%f\t%f\t%f\n',climateOS(i,1),climateOS(i,2),climateOS(i,3),climateOS(i,4),climateOS(i,5),climateOS(i,6),climateOS(i,7));
end
fclose(fid);

%% CropMix
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['CropMix_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileCropM=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Crop mix options for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Number of crop options %%');
fprintf( fid, '%s\n','1');
fprintf( fid, '%s\n','%% Specified planting calendar %%');
fprintf( fid, '%s\n','N');
fprintf( fid, '%s\n','%% Crop rotation filename %%');
fprintf( fid, '%s\n','N/A');
fprintf( fid, '%s\n','%% Information about each crop type %%');
fprintf( fid, '%s\n','%% CropType, CropFilename, IrrigationFilename %%');
fprintf( fid, '%s\n',['Maize, ' fileCrop ', IrrigationManagement.txt']);
fclose(fid);

%% Clock
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['Clock_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileClock=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Clock parameter inputs for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Simulation start time (yyyy-mm-dd) %%');
fprintf( fid, '%s%s%s\n','SimulationStartTime : ',dates_mvl{1},'');
fprintf( fid, '%s\n','%% Simulation end time (yyyy-mm-dd) %%');
fprintf( fid, '%s%s%s\n','SimulationEndTime : ',dates_mvl{2},'');
fprintf( fid, '%s\n','%% Simulate off-season ("N" or "Y") %%');
fprintf( fid, '%s\n','OffSeason : Y');
fclose(fid);

%% Soil
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['Soil_' num2str(j) '_' num2str(w.ProcessId) '_.txt';];
fileSoil=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Soil parameter inputs for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Soil profile filename %%');
fprintf( fid, '%s\n',['SoilProfile_' num2str(j) '_' num2str(w.ProcessId) '_.txt']);
fprintf( fid, '%s\n','%% Soil textural properties filename %%');
fprintf( fid, '%s\n',['SoilTexture_' num2str(j) '_' num2str(w.ProcessId) '_.txt']);
% fprintf( fid, '%s\n','N/A');
fprintf( fid, '%s\n','%% Soil hydraulic properties filename %%');
fprintf( fid, '%s\n','N/A');%SoilHydrology.txt
% fprintf( fid, '%s\n','SoilHydrology.txt');
fprintf( fid, '%s\n','%% Calculate soil hydraulic properties (0: No, 1: Yes) %%');
fprintf( fid, '%s\n','CalcSHP : 1'); %0
fprintf( fid, '%s\n','%% Total thickness of soil profile (m) %%');
fprintf( fid, '%s%f\n','zSoil : ',SoilInfo(1));
fprintf( fid, '%s\n','%% Total number of compartments %%');
fprintf( fid, '%s\n','nComp : 12');
fprintf( fid, '%s\n','%% Total number of layers %%');
fprintf( fid, '%s\n','nLayer : 1');
fprintf( fid, '%s\n','%% Thickness of soil surface skin evaporation layer (m) %%');
fprintf( fid, '%s\n','EvapZsurf : 0.04');
fprintf( fid, '%s\n','%% Minimum thickness of full soil surface evaporation layer (m) %%');
fprintf( fid, '%s\n','EvapZmin : 0.15');
fprintf( fid, '%s\n','%% Maximum thickness of full soil surface evaporation layer (m) %%');
fprintf( fid, '%s\n','EvapZmax : 0.30');
fprintf( fid, '%s\n','%% Maximum soil evaporation coefficient %%');
fprintf( fid, '%s\n','Kex : 1.1');
fprintf( fid, '%s\n','%% Shape factor describing reduction in soil evaporation %%');
fprintf( fid, '%s\n','fevap : 4');
fprintf( fid, '%s\n','%% Proportional value of Wrel at which soil evaporation layer expands %%');
fprintf( fid, '%s\n','fWrelExp : 0.4');
fprintf( fid, '%s\n','%% Maximum coefficient for soil evaporation reduction due to sheltering effect of withered canopy %%');
fprintf( fid, '%s\n','fwcc : 50');
fprintf( fid, '%s\n','%% Adjust default value for readily evaporable water (0: No, 1: Yes) %%');
fprintf( fid, '%s\n','AdjREW : 1'); %0
fprintf( fid, '%s\n','%% Readily evaporable water (mm) (only used if adjusting) %%');
fprintf( fid, '%s\n','REW : 13'); %8
fprintf( fid, '%s\n','%% Adjust curve number for antecedent moisture content (0:No, 1:Yes) %%');
fprintf( fid, '%s\n','AdjCN : 1');
fprintf( fid, '%s\n','%% Curve number %%');
fprintf( fid, '%s%f\n','CN : ',SoilInfo(3));
fprintf( fid, '%s\n','%% Thickness of soil surface (m) used to calculate water content to adjust curve number %%');
fprintf( fid, '%s\n','zCN : 0.3');
fprintf( fid, '%s\n','%% Thickness of soil surface (m) used to calculate water content for germination %%');
fprintf( fid, '%s\n','zGerm : 0.3');
fprintf( fid, '%s\n','%% Depth of restrictive soil layer (set to negative value if not present) %%');
fprintf( fid, '%s%f\n','zRes : ',SoilInfo(2)-0.001); %maak het net iets ondieper dan ST, anders fout code
fprintf( fid, '%s\n','%% Capillary rise shape factor %%');
fprintf( fid, '%s\n','fshape_cr : 16');
fclose(fid);

%% Soil Profile
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['SoilProfile_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileSoilP=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Soil profile discretisation for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% CompartmentNo	Thickness(m)	LayerNo %%');
for i =1:12
    fprintf( fid, '%i\t%f\t%i\n',i,SoilInfo(1)/12,1);
end
fclose(fid);

%% Soil Texture
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['SoilTexture_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileSoilT=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Soil layer texture for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% LayerNo   Thickness(m)   Sand(%)    Clay(%)  OrgMat(%)   DensityFactor %%');
fprintf( fid, '%i\t%f\t%i\t%i\t%i\t%i\n',1,SoilInfo(1),SoilInfo(4)*100,SoilInfo(5)*100,SoilInfo(6),1);
fclose(fid);

%% Soil Hydrology
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['SoilHydrology_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileSoilH=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Soil hydraulic properties for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% LayerNo  LayerThickness(m)  thS(m3/m3)    thFC(m3/m3)  thWP(m3/m3)   Ksat(mm/day) %%');
fprintf( fid, '%i\t%f\t%f\t%f\t%f\t%f\n',1,SoilInfo(1),0.52,0.44,0.23,120);
fclose(fid);

%% CO2
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['MaunaLoaCO2_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileCO2=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- CO2 concentration at Mauna Loa for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Year   CO2(ppm) %%');
fprintf( fid, '%i\t%f\n',2000,CO2atTime);
fprintf( fid, '%i\t%f\n',2001,CO2atTime);
fclose(fid);

%% Irrigation
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['IrrigationManagement_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileIrr=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Irrigation management parameters for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Irrigation time-series filename %%');
fprintf( fid, '%s\n','N/A');
fprintf( fid, '%s\n','%% Irrigation method (0: Rainfed; 1: Soil moisture based; 2: Fixed interval; 3: Specified time series; 4: Net calculation) %%');
fprintf( fid, '%s\n','IrrMethod : 0');
fprintf( fid, '%s\n','%% Irrigation interval in days (only used if IrrMethod = 2) %%');
fprintf( fid, '%s\n','IrrInterval : 3');
fprintf( fid, '%s\n','%% Soil moisture target in FAO56 growth stage one (% of total PAW below which irrigation is triggered) %%');
fprintf( fid, '%s\n','SMT1 : 70');
fprintf( fid, '%s\n','%% Soil moisture target in FAO56 growth stage two (% of total PAW below which irrigation is triggered) %%');
fprintf( fid, '%s\n','SMT2 : 70');
fprintf( fid, '%s\n','%% Soil moisture target in FAO56 growth stage three (% of total PAW below which irrigation is triggered) %%');
fprintf( fid, '%s\n','SMT3 : 70');
fprintf( fid, '%s\n','%% Soil moisture target in FAO56 growth stage four (% of total PAW below which irrigation is triggered) %%');
fprintf( fid, '%s\n','SMT4 : 70');
fprintf( fid, '%s\n','%% Maximum irrigation depth (mm) %%');
fprintf( fid, '%s\n','MaxIrr : 15');
fprintf( fid, '%s\n','%% Irrigation application efficiency (%) %%');
fprintf( fid, '%s\n','AppEff : 90');
fprintf( fid, '%s\n','%% Net irrigation threshold moisture level (% of total PAW that will be maintained) %%');
fprintf( fid, '%s\n','NetIrrSMT : 80.5');
fprintf( fid, '%s\n','%% Percentage of soil surface wetted by irrigation %%');
fprintf( fid, '%s\n','WetSurf : 100');
fclose(fid);

%% Irrigation schedule
% Niet nodig, wordt niet gelezen
%% Initial conditions
% path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
file=['InitialWaterContent_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileIni=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Initial soil water content for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Type of value (''Prop'' (i.e. WP/FC/SAT), ''Num'' (i.e. XXX m3/m3), ''Pct'' (i.e. % TAW)) %%');
fprintf( fid, '%s\n','Num');
fprintf( fid, '%s\n','%% Method (''Depth'': Inteprolate depth points; ''Layer'': Constant value for each soil layer) %%');
fprintf( fid, '%s\n','Layer');
fprintf( fid, '%s\n','%% Number of input points (NOTE: Must be at least one point per soil layer) %%');
fprintf( fid, '%s\n','1');
fprintf( fid, '%s\n','%% Input data points (Depth/Layer   Value) %%');
fprintf( fid, '%s\t%f\n','1',SWCtemp);
fclose(fid);

% % path='C:\Data_Maarten\software\AQ_OSmvl2\AquaCropOS_v50a\Input\';
% file=['InitialWaterContent_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
% fileIni=file;
% filename=[path file];
% fid = fopen( filename, 'wt' );
% fprintf( fid, '%s\n','%% ---------- Initial soil water content for AquaCropOS ---------- %%');
% fprintf( fid, '%s\n','%% Type of value (''Prop'' (i.e. WP/FC/SAT), ''Num'' (i.e. XXX m3/m3), ''Pct'' (i.e. % TAW)) %%');
% fprintf( fid, '%s\n','Num');
% fprintf( fid, '%s\n','%% Method (''Depth'': Inteprolate depth points; ''Layer'': Constant value for each soil layer) %%');
% fprintf( fid, '%s\n','Depth');
% fprintf( fid, '%s\n','%% Number of input points (NOTE: Must be at least one point per soil layer) %%');
% fprintf( fid, '%s\n','3');
% fprintf( fid, '%s\n','%% Input data points (Depth/Layer   Value) %%');
% fprintf( fid, '%s\t%f\n','0',SWCtemp);
% fprintf( fid, '%s\t%f\n','0.5',SWCtemp);
% fprintf( fid, '%s\t%f\n','0.55',0);
% fclose(fid);

%% Water table
file=['WaterTable_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileWat=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% --------- Groundwater table for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Water table present (''Y'' or ''N'') %%');
fprintf( fid, '%s\n','N');
fprintf( fid, '%s\n','%% Water table method (''Constant'' depth; ''Variable'' depth) %%');
fprintf( fid, '%s\n','Constant');
fprintf( fid, '%s\n','%% Date (dd/mm/yyyy)   Depth (m) %%');
fprintf( fid, '%s\n','01/09/2000  2');
fclose(fid);

%% Field Management
file=['FieldManagement_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
fileFieldM=file;
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ---------- Soil parameter inputs for AquaCropOS ---------- %%');
fprintf( fid, '%s\n','%% Soil surface covered by mulches (0: No; 1: Yes) %%');
fprintf( fid, '%s\n','Mulches : 0');
fprintf( fid, '%s\n','%% Area of soil surface covered by mulches during growing season (%) %%');
fprintf( fid, '%s\n','MulchPctGS : 50');
fprintf( fid, '%s\n','%% Area of soil surface covered by mulches outside growing season (%) %%');
fprintf( fid, '%s\n','MulchPctOS : 50');
fprintf( fid, '%s\n','%% Soil evaporation adjustment factor due to effect of mulches %%');
fprintf( fid, '%s\n','fMulch : 0.5');
fprintf( fid, '%s\n','%% Surface bunds present (0: No; 1: Yes) %%');
fprintf( fid, '%s\n','Bunds : 0');
fprintf( fid, '%s\n','%% Bund height (m) %%');
fprintf( fid, '%s\n','zBund : 0');
fprintf( fid, '%s\n','%% Initial water height in surface bunds (mm) %%');
fprintf( fid, '%s\n','BundWater : 0');
fclose(fid);

%% FileSetup
file=['FileSetup_' num2str(j) '_' num2str(w.ProcessId) '_.txt'];
filename=[path file];
fid = fopen( filename, 'wt' );
fprintf( fid, '%s\n','%% ----------- Input filenames for AquaCropOS ----------- %%');
fprintf( fid, '%s\n','%% Clock Filename %%');
fprintf( fid, '%s\n', fileClock);
fprintf( fid, '%s\n','%% Weather Data Filename %%');
fprintf( fid, '%s\n',fileClim);
fprintf( fid, '%s\n','%% C02 Concentration Time-Series Filename %%');
fprintf( fid, '%s\n',fileCO2);
fprintf( fid, '%s\n','%% Crop Parameters Filename %%');
fprintf( fid, '%s\n',fileCropM);
fprintf( fid, '%s\n','%% Soil Parameters Filename %%');
fprintf( fid, '%s\n',fileSoil);
fprintf( fid, '%s\n','%% Field Management Filename %%');
fprintf( fid, '%s\n',fileFieldM);
fprintf( fid, '%s\n','%% Initial Water Content Filename %%');
fprintf( fid, '%s\n',fileIni);
fprintf( fid, '%s\n','%% Groundwater Table Filename %%');
fprintf( fid, '%s\n',fileWat);
fprintf( fid, '%s\n','%% Output Filename %%');
fprintf( fid, '%s\n',['sample_' num2str(j)]);
fprintf( fid, '%s\n','%% Write daily outputs (''Y'' or ''N'') %%');
fprintf( fid, '%s\n','Y');
fclose(fid);
end




