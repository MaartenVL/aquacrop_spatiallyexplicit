function [BM,CY,Irriga,Infill,Drain,Runoff,ET,SWC,SWCfull,Evapo,Trans,StressW,StressT,GDcount,ImaxTot]=AquaCropOS_RUN(climateOS,dates_mvl,SoilInfo,CO2atTime,SWCtemp,AQpath,j,w,CropReadFile,CropReadFileTrees,LUval,Slope_input,climateOS_Pori,extraP4cell)
% ---------------------------------------------------------------------- %
% Tim Foster                                                             %
% June 2016                                                              %
%                                                                        %
% Function to run AquaCrop-OS v5.0a                                      %
%                                                                        %
% ---------------------------------------------------------------------- %

%% Declare global variables %%
global AOS_ClockStruct

%% Write text files
write_mvl_input(climateOS,dates_mvl,SoilInfo,CO2atTime,SWCtemp,AQpath,j,w,CropReadFile,CropReadFileTrees,LUval);
%% Run model %%
% Initialise simulation
AOS_Initialize(j,w);

% Perform single time-step (day)
while AOS_ClockStruct.ModelTermination == false
   AOS_PerformTimeStep(Slope_input,LUval,climateOS_Pori,extraP4cell);
end

% Finish simulation
AOS_Finish();
%% Output for 1 year
global AOS_InitialiseStruct
%% write yearly export variables
BM=max(AOS_InitialiseStruct.Outputs.CropGrowth(:,11));
CY=max(AOS_InitialiseStruct.Outputs.CropGrowth(:,15));
Irriga=AOS_InitialiseStruct.Outputs.WaterFluxes(:,9);
Infill=AOS_InitialiseStruct.Outputs.WaterFluxes(:,10);
Drain=AOS_InitialiseStruct.Outputs.WaterFluxes(:,12);
Runoff=AOS_InitialiseStruct.Outputs.WaterFluxes(:,11);
ET=AOS_InitialiseStruct.Outputs.WaterFluxes(:,15)+AOS_InitialiseStruct.Outputs.WaterFluxes(:,17);
SWC=mean(AOS_InitialiseStruct.Outputs.WaterContents(end-1,6:17));
%SWCfull=AOS_InitialiseStruct.Outputs.WaterContents(:,6:17);
SWCfull=mean(mean(AOS_InitialiseStruct.Outputs.WaterContents(:,6:17)));
Evapo=AOS_InitialiseStruct.Outputs.WaterFluxes(:,15);
Trans=AOS_InitialiseStruct.Outputs.WaterFluxes(:,17);
%StressW=AOS_InitialiseStruct.Outputs.Stresses.Water;% MEMORY
%StressT=AOS_InitialiseStruct.Outputs.Stresses.Temp;% MEMORY
%StressW=[StressW;ones(size(StressW(1,:)))]; % onderste rij missing anders (364 rijen ipv 365)
%StressT=[StressT;ones(size(StressT(1,:)))];
% getting % stress
%StressW=(StressW-ones(size(StressW)))*-1;
%StressT=(StressT-ones(size(StressT)))*-1;
%GDcount=AOS_ClockStruct.GDcount;
StressW=0;
StressT=0;
GDcount=0;
%Imax
Imax=AOS_InitialiseStruct.Outputs.Imax;
ImaxTot=sum(Imax);
%% Remove Input and Output files (with parfor index j)
delete([AQpath 'Input/*' num2str(j) '_' num2str(w.ProcessId) '_.txt']);
delete([AQpath 'Output/*' num2str(j) '_' num2str(w.ProcessId) '_.txt']);
end