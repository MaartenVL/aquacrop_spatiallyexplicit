function main_rout4os_WF(time)

%% convert time to number
time=str2num(time);

%% *** Initialising - General - add paths ***
tic
tak=0; % counter

%Trees4AQ
tree4aq=0; % to model effect trees differently - mostly switched off

% Rainfall data source
check_Rainfall=1; % 0 = Hans Renssen - 1 = WeaGETS

%% Path's
datapath='/scratch/leuven/318/vsc31850/AQ_on_HPC_parallel_v7/AQ_OSmvl2/input_mvl_files/';
AQpath='/scratch/leuven/318/vsc31850/AQ_on_HPC_parallel_v7/AQ_OSmvl2/AquaCropOS_v50a/';

% write Filelocations text file
write_fileLoc(AQpath,time);


%% *** Read Input data ***

%DTM
file = sprintf('v4_5_5_carved_idri');
fname= [datapath file];
[dtm,dtmrow,dtmcol,dtmhres,dtmvres] = Read_Idrisi_dtm(fname);

% Slope
file='slope_percent.tif';
fname= [datapath file];
slopemap=imread(fname)/100;
slopemap=double(slopemap);

% climate data -  Hans Renssen Temp climate data
file='T_HR.txt'
filename = [datapath file];
delimiter = '\t';
formatSpec = '%*s%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
T_HR = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;
file='corrfactTHR_no06_corr.mat';
fname=[datapath file];
corrfact=load(fname);
T_HR=T_HR.*repmat(corrfact.corr_fact,4000,1);

% Temp data omzetten Tmin & Tmax
%based on form. derived from NewLocClim (see C:\Data_Maarten\Analysis\Soil
%Productivity\gravgaz\Gravgaz input parameters_nieuw)
Tmin_in=(0.8533.*T_HR)-3.4119;
Tmax_in=(1.1463.*T_HR)+3.8445;

% inlezen Hans Renssen Precip climate data
if check_Rainfall == 0
    file='P_HR_AQ.txt';
    filename = [datapath file];
    delimiter = '\t';
    formatSpec = '%*s%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    P_HR = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans;
    
    % Randomly generated daily data via Poesen IDF curves
    file='P_HR_randPoes.mat';
    fname=[datapath file];
    load(fname);
    P_HR_rand=P_Poes;
    
    % WeaGETS - day 1 = 1 january xxxx !
else
    file='gP_P_HR.mat';
    fname=[datapath file];
    load(fname);
    P_HR_rand=gP;
end

% SWC multiple regress
swc_param=load([datapath 'swc_param.mat']);
swc_param=swc_param.b;
swc_paramT=load([datapath 'swc_paramT.mat']);
swc_paramT=swc_paramT.bTree;

%% *** selection sub grid ***

xx1=37;% full ctch
xx2=264; % full ctch
yy1=93;% full ctch
yy2=344;% full ctch
nx=xx2-xx1;
ny=yy2-yy1;


%% *** TESTMODI ***
testmodeLU=0;
testmodeLUshort=0;
testmodeP=0;
testmodeST=0;
testmodeSTshort=0;%1= testmode, 2=1e 20 jaar elke jaar run
testmodeDTM=0;

%testveld N-Marsh
if testmodeDTM==1
    xx1=147;
    xx2=160;
    yy1=131;
    yy2=153;
    nx=xx2-xx1;
    ny=yy2-yy1;
end

dtm2=dtm(yy1-1:yy2+1,xx1-1:xx2+1); %om de randpixels van dtm een juiste afvoer te geven
dtm=dtm(yy1:yy2,xx1:xx2);


%testP
file='Ptest.mat';
fname=[datapath file];
Pt=load(fname, 'Pday365')
Ptest=Pt.Pday365;
%testLU
testLUdeg=ones(374,349);
testLUundeg=ones(374,349)+1;
%STtest
testST=ones(nx,ny)+900;

%% *** Lithology & other GIS data ***
%lithology
file='lithology.tif';
fname=[datapath file];
lithology=imread(fname);

%marsh area
file='cmap1.tif';
fname=[datapath file];
landuse4mask=imread(fname);
lithologymask=lithology;
lithologymask(lithology>0)=1;
marshmask=lithology;
marshmask(:,:)=0;
marshmask(lithology > 0 & landuse4mask == 0)=1;

%small Marsh
file='marshsmallr1.tif';
fname=[datapath file];
marshS=imread(fname);
marshS=double(marshS);
marshS=int16(marshS);
marshS(marshS==0)=10;
marshS(marshS==255)=0;
marshSuper=marshS+marshmask;
marshSuper(marshSuper<11)=0;
marshSuper(marshSuper==11)=1;

%sedexport
file='sedexport.mat';
fname=[datapath file];
sedexport=load(fname);
sedexport=sedexport.odtotalponddeprivexpor;
sedexport=sedexport.*10; % since 10 years for each time step that accumulate this sedexport to marsh
sedexport=cumsum(sedexport); % sediments accumulate, do not dissapear.
%size grid
nrow=numel(dtm(:,1));
ncol=numel(dtm(1,:));

%% *** Initialization ***
SWCtemp(:,:)=0;

%% *** 0.Time Period ***
[dates_mvl,dates4climOS]=get_time_period;

%% *** 2.Soil & Crop Load ***
ParamStruct=readCrop4os(AQpath);

%3.Irrigation
%IrrFile=write_irr(ST);

%4.Management
%ManaFile=write_mana();

%9.DayFile
% [FirstDaySim,LastDaySim,FirstDayCrop,LastDayCrop]=write_dayfile_yearly(FirstDaySimRead,LastDaySimRead,FirstDayCropRead,LastDayCropRead,FirstMonSimRead,LastMonSimRead,FirstMonCropRead,LastMonCropRead)

%% TIME
maxtime=400; %400
dtime=2; %20

%% SINGLE TIME
singletime=0;
if singletime==1
    maxtime=400;
    dtime=1;
end

%Output initialising
OutputSummarytext={'ST (mm)','CY (t/ha/yr)','BioMass (kg/ha/yr)','Infill (mm)','Drain (mm)','Runoff at pixel (mm)','Precip (mm)','Precip + Incomming Runoff (mm)','Evapo (mm)','Trans (mm)','SWC at end sim(mm)','Export (mm)','Irri (mm)','CO2 (ppm)'};
OutputSummary=zeros(1,numel(OutputSummarytext));

%% *** time ***
%for time=1:1 % SINGLE
     %time_actual=time_array(time); %uit single
%time_actual=400; %SINGLE
%LU_timer=5; %SINGLE
    if time < 110
        LU_timer=1;
    elseif time < 130
        LU_timer=2;
    elseif time < 150
        LU_timer=3;
    elseif time < 210
        LU_timer =4;
    elseif time < 310
        LU_timer=5;
    elseif time < 380
        LU_timer=6;
    elseif time < 390
        LU_timer=5;
    else
        LU_timer=5;
    end
    if testmodeLUshort == 1;
        if time < 6
            LU_timer=1;
        elseif time < 7
            LU_timer=2;
        elseif time < 8
            LU_timer=3;
        elseif time < 11
            LU_timer =4;
        elseif time < 16
            LU_timer=5;
        elseif time < 19
            LU_timer=6;
        elseif time < 20
            LU_timer=5;
        end
    end
    %LandUse
    file=sprintf('%s%d%s','cmap',LU_timer,'.tif');
    fname= [datapath file];
    landuse=imread(fname);
    landuse(landuse==20000)=1; % degraded
    landuse(landuse==10000)=2; % undegraded
    landuse=landuse+marshmask;
    if testmodeLU ==1
        landuse=testLUdeg;
    end
    landuse(marshSuper==1)=0;
    
    % Routing algorithm
    [runoffcalc,rowmatrix,colmatrix] = TopoRouting(dtm,dtm2,nrow,ncol,AQpath); % see TopoRouting Toolbox: Schwanghart, W. and Kuhn, N. J. (2010), Environmental Modelling & Software TopoToolbox
    [waarde,locatie]=sort(runoffcalc(:),'ascend');
    [rowY,colX]=ind2sub(size(runoffcalc),locatie);
    

    nx=xx2-xx1+1; %the extent of x
    ny=yy2-yy1+1; %the extent of y
    outerPix=2*ny+(2*(nx-2));
    rowY(1:outerPix)=[];
    colX(1:outerPix)=[];
    
    %1.Climate File
    [Tmin,Tmax,T,P_HR_randout,Pprev]=read_T_P_yearly_mvl(time,dtime,T_HR,P_HR_rand,Tmax_in,Tmin_in); %P should be daily so that runoff is influenced by water content soil (or CN influenced by SWC)
    [EToday,Rn,Tday,TmaxDay,TminDay]=write_ETo_mvl_v2(Tmin,Tmax,T); % ETo and Tday interpolation
    CO2atTime=write_CO2_mvl(time,dtime,datapath);

    PMat=P_HR_randout;
    climateOS(:,1)=TminDay;
    climateOS(:,2)=TmaxDay;
    climateOS(:,3)=P_HR_randout;
    climateOS(:,4)=EToday;
    climateOS=[dates4climOS climateOS];
    
    %AQ - 1 month (time in main, space in AQ batcher)
    [Runtemp,CYtemp,BMtemp,Irrtemp,Itemp,STmap,Draintemp,SWCtemp,extraP4cellSum,Pouttemp,InTree,extraP4cell2D,Evapotemp,Transtemp,StressVars,SWCfull_temp,GDcount_temp,CN_temp,extraP4cell_days_temp,check_rout_P_extraP,SWCtemp_INPUT,Imaxtemp]=AQ_rout_mvl(time, xx1, xx2, yy1, yy2,nrow,ncol,dates_mvl,climateOS,landuse,lithology,sedexport,marshmask,rowY,colX,rowmatrix,colmatrix,outerPix,SWCtemp,testmodeSTshort,tree4aq,CO2atTime,ParamStruct,datapath,AQpath,runoffcalc,slopemap,swc_param,Pprev,swc_paramT);

    %% Handling results
    STmapStore=STmap;
    CY=CYtemp;
    BM=BMtemp;
    Irr=Irrtemp;
    I=Itemp;
    Drain=Draintemp;
    RunStore=Runtemp;
    extraP4cell2Dstore=extraP4cell2D;
    EvapoStore=Evapotemp;
    TransStore=Transtemp;
    EvapoStoreYear=sum(EvapoStore,3);
    TransStoreYear=sum(TransStore,3);
    Run2D=sum(RunStore,3);
    CN_store=CN_temp;
    check_rout_P_extraP_store=check_rout_P_extraP;
    check_rout_P_extraPsum_store=sum(sum(check_rout_P_extraP));
    SWCinput=SWCtemp_INPUT;
    Imax=Imaxtemp;
      
    RunMat=squeeze(sum(sum(Runtemp(:,:,:)))); %runoff lezen we elke dag
    IMat=sum(sum(Itemp(:,:,:))); %infil enkel totaal op einde jaar
    DrainMat=sum(sum(Draintemp(:,:,:))); % Drain enkel totaal op einde jaar
    %extraP4cellstore(time,:)=sum(extraP4cell);
    %extraP4cellstoreall(time,:)=extraP4cell;
    extraP4cellSumYear=sum(extraP4cellSum);
    extraP4cellSumStore=extraP4cellSum;
    extraP4cellMat=squeeze(sum(sum(extraP4cell_days_temp)));
    extraP4cell_days_tempStore=extraP4cell_days_temp;
    extraP4cell_days_temp_Total=sum(sum(sum(extraP4cell_days_temp)));
    PoutMat_temp=sum(sum(Pouttemp(:,:,:))); % P + extra P
    PoutMat=PoutMat_temp;
    %     waterexport(time)=(sum(P)+sum(RunMat(time,:)))-(sum(sum(extraP4cellstore{time})));
    %     waterexport(time)=(sum(P)+sum(RunMat(time,:)))-(extraP4cellstore(time));
     waterexport_temp=sum(RunMat)-sum(sum(sum(extraP4cell_days_temp)));
    waterexport=waterexport_temp;
    waterexport_day=RunMat-squeeze(sum(sum(extraP4cell_days_temp)));
    SWCstore=SWCtemp;
    InMat=sum(sum(InTree(:,:,:)));
    

     SWCVars2Dstore=SWCfull_temp;
     

landusedomain=landuse(yy1+1:yy2-1,xx1+1:xx2-1);
landusedomain(landusedomain > 0)=1;
landsum=sum(sum(landusedomain));
PMatSum=sum(PMat,2)*landsum; % Met grote delen buiten ctc gaat dit fout zijn de *nx-1 ny-1 methode
checkPandRun=PoutMat-PMatSum-extraP4cell_days_temp_Total';
checkP_I_R=PoutMat-sum(RunMat)-IMat;
checkwatbal=waterexport'+extraP4cell_days_temp_Total'+IMat-PoutMat;
check_Store{1}=checkPandRun;
check_Store{2}=checkP_I_R;
check_Store{3}=checkwatbal;
check_Store{4}=check_rout_P_extraPsum_store;


RunTimeStored=toc



% clock
c=clock;


c_cell=num2cell(c);
Tval=cellfun(@num2str,c_cell,'UniformOutput',0)

% Create directory for storage
% -----------------------------
dirname=[AQpath 'Output/' 'worker' num2str(time) '_' Tval{1} '_' Tval{2} '_' Tval{3} '_' Tval{4} '_' Tval{5} '_' Tval{6} '_Final'];
mkdir(dirname);

% Store to Disk
% --------------
save([dirname '/' 'CY.mat'],'CY');
%save([dirname '/' 'STmapStore.mat'],'STmapStore');
save([dirname '/' 'BM.mat'],'BM');
save([dirname '/' 'I.mat'],'I');
save([dirname '/' 'Drain.mat'],'Drain');
save([dirname '/' 'EvapoStoreYear.mat'],'EvapoStoreYear');
save([dirname '/' 'TransStoreYear.mat'],'TransStoreYear');
save([dirname '/' 'SWCstore.mat'],'SWCstore'); %2D - end
save([dirname '/' 'Run2D.mat'],'Run2D');
%save([dirname '/' 'P_HR_rand.mat'],'P_HR_rand');
save([dirname '/' 'PoutMat.mat'],'PoutMat');
save([dirname '/' 'waterexport.mat'],'waterexport');
save([dirname '/' 'extraP4cell2Dstore.mat'],'extraP4cell2Dstore');
%save([dirname '/' 'StressVarsStore' num2str(time) '.mat'],'StressVarsStore');
save([dirname '/' 'SWCVars2Dstore.mat'],'SWCVars2Dstore'); %2D - mean
save([dirname '/' 'RunMat.mat'],'RunMat');
%save([dirname '/' 'landusestore.mat'],'landusestore');
%save([dirname '/' 'runoffcalc.mat'],'runoffcalc');
save([dirname '/' 'CN_store.mat'],'CN_store');
save([dirname '/' 'extraP4cellMat.mat'],'extraP4cellMat');
save([dirname '/' 'check_rout_P_extraP_store.mat'],'check_rout_P_extraP_store');
save([dirname '/' 'RunTimeStored.mat'],'RunTimeStored');
save([dirname '/' 'check_Store.mat'],'check_Store');
%save([dirname '/' 'OutFile.mat'],'OutFile');
save([dirname '/' 'SWCinput.mat'],'SWCinput'); %2D - start
save([dirname '/' 'Imax.mat'],'Imax');

