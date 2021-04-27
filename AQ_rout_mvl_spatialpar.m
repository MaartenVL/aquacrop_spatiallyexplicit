function [Runtemp,CYtemp,BMtemp,Irrtemp,Itemp,STmap,Draintemp,SWCtemp,extraP4cellSum,Pouttemp,InTree,extraP4cell2D,Evapotemp,Transtemp,StressVars,SWCVars,GDcount_temp,CN_temp,extraP4cell_days_temp,check_rout_P_extraP] = AQ_rout_mvl(time, xx1, xx2, yy1, yy2,nrow,ncol,dates_mvl,climateOS,landuse,lithology,sedexport,marshmask,rowY,colX,rowmatrix,colmatrix,outerPix,SWCtemp,testmodeSTshort,tree4aq,CO2atTime,ParamStruct,datapath,AQpath,runoffcalc,time_actual)
% reading crop file
% CropReadFile=load([datapath 'CropRead.mat']); % old, not calibrated
CropReadFile=load([datapath 'cropfile_cal4paleo.mat']); % new, calibrated values of winter wheat
CropReadFileTrees=load([datapath 'cropfile_cal4paleo4trees.mat']); % new, calibrated values of winter wheat
CropReadFile.CropRead=CropReadFile.CropwheatCal4Paleo;

%reading mask file
% [mask, R] = imread('C:\Data_Maarten\Analysis\watem sedem\MATLAB for WS\STmap\V4_5_5_PRC_ML_tif.tif');
% mask=ones(4,4); %to do (actually read land use) when it's an agricult. field, we want to run AQ, later this comes back as mask = 1? ok run AQ. if zero we don't need to know the CY

statusMrow=[1,1,1;0,0,0;-1,-1,-1];
statusMcol=[1,0,-1;1,0,-1;1,0,-1];
% extraP=zeros(yy2-yy1+1,xx2-xx1+1,334);
extraPtemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
extraP4cellSum=0;
% Runtemp=zeros(yy2-yy1+1,xx2-xx1+1,334);
Runtemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
Evapotemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
Transtemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
Irrtemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
StressW_Aer_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);
StressW_Exp_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);
StressW_Pol_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);
StressW_Sen_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);
StressW_Sto_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);
StressT_Bio_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);
StressT_PolC_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);
StressT_PolH_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_1=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_2=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_3=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_4=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_5=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_6=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_7=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_8=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_9=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_10=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_11=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp_12=zeros(yy2-yy1+1,xx2-xx1+1,365);
extraP4cell_days_temp=zeros(yy2-yy1+1,xx2-xx1+1,365);

%%
CYtemp=zeros(yy2-yy1+1,xx2-xx1+1);
Itemp=zeros(yy2-yy1+1,xx2-xx1+1);
Draintemp=zeros(yy2-yy1+1,xx2-xx1+1);
BMtemp=zeros(yy2-yy1+1,xx2-xx1+1);
STmap=zeros(yy2-yy1+1,xx2-xx1+1);
Pouttemp=zeros(yy2-yy1+1,xx2-xx1+1);
InTree=zeros(yy2-yy1+1,xx2-xx1+1);
extraP4cell2D=zeros(yy2-yy1+1,xx2-xx1+1);
GDcount_temp=zeros(yy2-yy1+1,xx2-xx1+1);
CN_temp=zeros(yy2-yy1+1,xx2-xx1+1);

if time == 1
    SWCtemp=zeros(yy2-yy1+1,xx2-xx1+1);
end
%% Preparation parfor loop
% [buurrow,buurcol,rowmatrix,colmatrix,runoffcalc] = routing(dtm,landuse,nrow,ncol,dtm2);
runoffcalcmasked=runoffcalc;
runoffcalcmasked([1 end],:)=0;
runoffcalcmasked(:,[1 end])=0;
[waarde,locatie]=sort(runoffcalcmasked(:),'ascend');
[rowY,colX]=ind2sub(size(runoffcalcmasked),locatie);

%Removing border pixels with 0 cells receiving (if not removed, they will
%throw an index error
nx=xx2-xx1+1; %the extent of x
ny=yy2-yy1+1; %the extent of y
outerPix=2*ny+(2*(nx-2));
rowY(1:outerPix)=[]; % omdat we geen 0 runoff hebben, outerPix=/= de randen
colX(1:outerPix)=[];
runoffcalcmasked=runoffcalcmasked(2:end-1,2:end-1);

%% initialize output matrices
% calcresult_ini=zeros(size(runoffcalcmasked));
uniquerout=unique(runoffcalcmasked);
num_it=numel(uniquerout);
for i=1:num_it
    it_num=uniquerout(i);
    it_num_amount=numel(runoffcalcmasked(runoffcalcmasked==it_num));
    if i ==1
        colX_select_pool=colX;
        rowY_select_pool=rowY;
    end
    colX_select=colX_select_pool(1:it_num_amount);
    rowY_select=rowY_select_pool(1:it_num_amount);

    % Initialize SLICE variables
    STmap_SLICE=zeros(1,it_num_amount);
    Runtemp_SLICE=zeros(365,it_num_amount);
    extraP4cell2D_SLICE=zeros(1,it_num_amount);
    CYtemp_SLICE=zeros(1,it_num_amount); % reden dat deze vars 6x6 zijn, is omdat de Linkse top en zijde 0 krijgen omdat 1e cel is op 2,2, en we later nog extra rij en col bij voegen rechts zijde en onder
    %extraP_SLICE=zeros(365,it_num_amount);
    Irrtemp_SLICE=zeros(365,it_num_amount);
    Itemp_SLICE=zeros(1,it_num_amount);
    Draintemp_SLICE=zeros(1,it_num_amount);
    BMtemp_SLICE=zeros(1,it_num_amount);
    SWCtemp_SLICE=zeros(1,it_num_amount);
    Pouttemp_SLICE=zeros(1,it_num_amount);
    InTree_SLICE=zeros(1,it_num_amount);
    Evapotemp_SLICE=zeros(365,it_num_amount);
    Transtemp_SLICE=zeros(365,it_num_amount);
    StressW_Aer_SLICE=zeros(365,it_num_amount);
    StressW_Exp_SLICE=zeros(365,it_num_amount);
    StressW_Pol_SLICE=zeros(365,it_num_amount);
    StressW_Sen_SLICE=zeros(365,it_num_amount);
    StressW_Sto_SLICE=zeros(365,it_num_amount);
    StressT_Bio_SLICE=zeros(365,it_num_amount);
    StressT_PolC_SLICE=zeros(365,it_num_amount);
    StressT_PolH_SLICE=zeros(365,it_num_amount);
    SWCfull_SLICE_1=zeros(365,it_num_amount);
    SWCfull_SLICE_2=zeros(365,it_num_amount);
    SWCfull_SLICE_3=zeros(365,it_num_amount);
    SWCfull_SLICE_4=zeros(365,it_num_amount);
    SWCfull_SLICE_5=zeros(365,it_num_amount);
    SWCfull_SLICE_6=zeros(365,it_num_amount);
    SWCfull_SLICE_7=zeros(365,it_num_amount);
    SWCfull_SLICE_8=zeros(365,it_num_amount);
    SWCfull_SLICE_9=zeros(365,it_num_amount);
    SWCfull_SLICE_10=zeros(365,it_num_amount);
    SWCfull_SLICE_11=zeros(365,it_num_amount);
    SWCfull_SLICE_12=zeros(365,it_num_amount);
    GDcount_SLICE=zeros(1,it_num_amount);
    CN_SLICE=zeros(1,it_num_amount);
    extraP4cell_days_SLICE=zeros(365,it_num_amount);

  %parpool(100)
    parfor j=1:it_num_amount
         %j
        % w.ProcessId=505;
        w=getCurrentWorker;
        %disp(j);
        % for i=1:nrow*ncol
        %     tel=tel+1
        %         ii=colX(i);
        %         jj=rowY(i);
        ii=colX_select(j);
        jj=rowY_select(j);
        %disp(['Slave ' num2str(j) ' : ' 'just before if statement']);
        if landuse(jj+yy1-1,ii+xx1-1) == 1 || landuse(jj+yy1-1,ii+xx1-1) == 2; % yy1 en xx1 zorgen ervoor dat je de juiste LU pixels selecteert
            LUval=landuse(jj+yy1-1,ii+xx1-1);% LUval for crop decision
             %disp(['Slave ' num2str(j) ' : ' 'entered if statement']);
            [ST]=get_soil_rout(time_actual,ii,jj,xx1,yy1,sedexport,marshmask,testmodeSTshort,datapath); %TO DO: fix time problem.% TIME SOIL =/= TIME RAIN. Soil in 10yearly, climate = yearly % <<--- fixed the last one
            STmap_SLICE(j)=ST;
            CN=write_CN(landuse,ST,lithology,jj,ii,xx1,yy1,datapath);
            CN_SLICE(j)=CN;

            controlMrow=rowmatrix(rowY_select(j)-1:rowY_select(j)+1,colX_select(j)-1:colX_select(j)+1); % cel die je beschouwd, die zijn buren, naar waar geven ze hun afvoer? Y richting
            controlMcol=colmatrix(rowY_select(j)-1:rowY_select(j)+1,colX_select(j)-1:colX_select(j)+1); % cel die je beschouwd, die zijn buren, naar waar geven ze hun afvoer? X richting

            diffMrow=statusMrow-controlMrow;
            diffMcol=statusMcol-controlMcol;

            binaryM=zeros(size(statusMrow));
            binaryM(diffMrow == 0 & diffMcol ==0)=1;

            extraPM=extraPtemp(rowY_select(j)-1:rowY_select(j)+1,colX_select(j)-1:colX_select(j)+1,:); %what if this pixel is outside ctc? it shouldnt matter since its runoff will be zero and won't affect the sum of extraP (see RunTemp, outer pixels are empty)

            binaryM=repmat(binaryM,[1,1,size(extraPM,3)]);
            extraP4cell=binaryM.*extraPM;
            extraP4cell=sum(sum(extraP4cell)); % 1e en 2e sum = x en y (2D), 3e dim blijft staan
            extraP4cell=extraP4cell(:)'; %dagelijkse extra neerslag voor bepaald jaar, voor bepaalde cell (via runoff van omliggende pixels)
            extraP4cellSum=extraP4cellSum+extraP4cell; %om op te tellen bij P voor de waterbalans => werkt dit in parfor?!

            % Climate will be different for each parfor j iteration, and
            % cant be a variable w/o indexing, since it will have a
            % different outcome depending on the iteration (climateOS=1 but
            % also climateOS=2)
            [climateOS_INPUT,Pout]=write_rain_rout365_mvl(climateOS,time,extraP4cell);
            %         ClimateFile=write_climate(TempFile,EToFile,RainFile,CO2File,time);

            %5.Soil File
            [ST,ST2]=write_soil_rout_mvl(ST);
            %read the texture map
            texture=get_texture4os(lithology,jj,ii,xx1,yy1);
            ParamStruct_INPUT=readSoil4os(ParamStruct,ST,ST2,CN,texture); % same reasoning as climateOS, altough we don't need ParamStruct since necessary data is written to SoilInfo
            SoilInfo=[ST ST2 CN texture];
            %SoilInfo(2,:)=texture;
            %6.Groundwater
            %GWFile=write_GW();

            %7.Initial Conditions

            if  time == 1; % enkel in de 1e run geen inicon file
                SWCtemp_INPUT=0.30;%voor SCL is dit de start FC=0.44
                %             SWCtemp=write_IniCon_rout_mvl(ST,ii,jj,SWCtemp); %SWCtemp na model run wordt doorgegeven naar volgend jaar, resultaat hiervan lezen we hier
            else
            %disp('time=/=1');
                SWCtemp_INPUT=SWCtemp(jj,ii);
            end
            %8.Off-seson conditions

            %         if tree4aq==1
            %             if landuse(jj+yy1-1,ii+xx1-1) == 1
            %                 if  time ~= 1;
            %                     write_projectfile_month_rout(FirstDaySim,LastDaySim,FirstDayCrop,LastDayCrop,ClimateFile,TempFile,EToFile,RainFile,CO2File,CropFile,SoilFile,IniConFile,ii,jj,time);
            %                 else
            %                     write_projectfile_rout(FirstDaySim,LastDaySim,FirstDayCrop,LastDayCrop,ClimateFile,TempFile,EToFile,RainFile,CO2File,CropFile,SoilFile,ii,jj,time);
            %                 end
            %             else
            %                 if time == 1
            %                     SWCtree=0;
            %                 else
            %                     SWCtree=SWCtemp(jj,ii);
            %                 end
            %                 [Rtree,Infiltree,Wlosstree,EWtree,InTreetic]=tree(EToday,P4tree,CN,ST,Rn,SWCtree,ii,jj,time);
            %                 CYtemp(jj,ii)=0;
            %                 Runtemp(jj,ii,:)=Rtree; % every day is saved!!!
            %                 Itemp(jj,ii)=Infiltree;
            %                 Draintemp(jj,ii)=Wlosstree;
            %                 BMtemp(jj,ii)=0;
            %                 SWCtemp(jj,ii)=EWtree;
            %                 Pouttemp(jj,ii)=Pout;
            %                 InTree(jj,ii)=InTreetic;
            %                 WatBal(jj,ii)=WatBaltic;
            %             end
            %
            %         else
            %10.Project file
            %             if  time ~= 1;
            %                 write_projectfile_month_rout(FirstDaySim,LastDaySim,FirstDayCrop,LastDayCrop,ClimateFile,TempFile,EToFile,RainFile,CO2File,CropFile,SoilFile,IniConFile,ii,jj,time);
            %             else
            %                 write_projectfile_rout(FirstDaySim,LastDaySim,FirstDayCrop,LastDayCrop,ClimateFile,TempFile,EToFile,RainFile,CO2File,CropFile,SoilFile,ii,jj,time);
            %             end
            %         end

            %Run AQ plugin
            %         [BM,CY,Irriga,Infill,Drain,Runoff,ET,SWC,Evapo,Trans] = AquaCropOS_RUN(ParamStruct,climateOS,SWCtemp,dates_mvl,CO2atTime);
            %disp(['Slave ' num2str(j) ' : ' '"I''m Goin In"']);
            [BM,CY,Irriga,Infill,Drain,Runoff,ET,SWC,SWCfull,Evapo,Trans,StressW,StressT,GDcount]=AquaCropOS_RUN(climateOS_INPUT,dates_mvl,SoilInfo,CO2atTime,SWCtemp_INPUT,AQpath,j,w,CropReadFile,CropReadFileTrees,LUval);
            %disp(['Slave ' num2str(j) ' : ' '"I''m Out"']);
            %         system('C:\Data_Maarten\software\AquaCropBatch5\ACsaV50.exe'); %geeft ans = 0 als het uitgevoerd is!!

            % if time == 1 %forgot why I put this here?
            %     clearvars I
            % end
            % if AQflag ==3 %checken! niet enkel AQflag 3
            %     for time = 1:400:4000 %Model run time

            %         flag=0;                                              %\
            %         flag=read_mask_month(mask,flag,time,i,j,xx1,yy1);    %| old, %before we explicetly took into account landuse
            %         if flag == 0;                                        %/
            %         if landuse(i+xx1,j+yy1) == 1 || landuse(i+xx1,j+yy1) == 2; %only when we are in the catchment I want to read the output, since we didnt calculate pixels OUTSIDE the catchment
            %visualize CY,R
            %map(jj,ii)=1;
            AQflag=3;%dit moet weg! geen nood aan AQflag in regel hieronder. haal dat uit de functie, voorlopig nog laten staan. er zijn 2 apparte leesfiles voor output, AQflag =/= nodig
            %         if tree4aq==1
            %             if landuse(jj+yy1-1,ii+xx1-1) == 1
            %                 [CYtic,Runtic,Itic,BMtic,Draintic,SWC]=read_output_rout(ii,jj,time,AQflag);
            %                 CYtemp(jj,ii)=CYtic; % reden dat deze vars 6x6 zijn, is omdat de Linkse top en zijde 0 krijgen omdat 1e cel is op 2,2, en we later nog extra rij en col bij voegen rechts zijde en onder
            %                 Runtemp(jj,ii,:)=Runtic;
            %                 Itemp(jj,ii)=Itic;
            %                 Draintemp(jj,ii)=Draintic;
            %                 BMtemp(jj,ii)=BMtic;
            %                 SWCtemp(jj,ii)=SWC;
            %                 Pouttemp(jj,ii)=Pout;
            %                 InTree(jj,ii)=0;
            %             end
            %         else
            %             [CYtic,Runtic,Itic,BMtic,Draintic,SWC]=read_output_rout(ii,jj,time,AQflag);
            CYtemp_SLICE(j)=CY; % reden dat deze vars 6x6 zijn, is omdat de Linkse top en zijde 0 krijgen omdat 1e cel is op 2,2, en we later nog extra rij en col bij voegen rechts zijde en onder
            Runtemp_SLICE(:,j)=Runoff;
            Irrtemp_SLICE(:,j)=Irriga;
            Itemp_SLICE(j)=sum(Infill);
            Draintemp_SLICE(j)=sum(Drain);
            BMtemp_SLICE(j)=BM;
            SWCtemp_SLICE(j)=SWC;
            Pouttemp_SLICE(j)=sum(Pout);
            InTree_SLICE(j)=0;
            extraP4cell2D_SLICE(j)=sum(extraP4cell);
            extraP4cell_days_SLICE(:,j)=extraP4cell;
            Evapotemp_SLICE(:,j)=Evapo;
            Transtemp_SLICE(:,j)=Trans;
            StressW_Aer_SLICE(:,j)=StressW(:,1);
            StressW_Exp_SLICE(:,j)=StressW(:,2);
            StressW_Pol_SLICE(:,j)=StressW(:,3);
            StressW_Sen_SLICE(:,j)=StressW(:,4);
            StressW_Sto_SLICE(:,j)=StressW(:,5);
            StressT_Bio_SLICE(:,j)=StressT(:,1);
            StressT_PolC_SLICE(:,j)=StressT(:,2);
            StressT_PolH_SLICE(:,j)=StressT(:,3);
            SWCfull_SLICE_1(:,j)=SWCfull(:,1);
            SWCfull_SLICE_2(:,j)=SWCfull(:,2);
            SWCfull_SLICE_3(:,j)=SWCfull(:,3);
            SWCfull_SLICE_4(:,j)=SWCfull(:,4);
            SWCfull_SLICE_5(:,j)=SWCfull(:,5);
            SWCfull_SLICE_6(:,j)=SWCfull(:,6);
            SWCfull_SLICE_7(:,j)=SWCfull(:,7);
            SWCfull_SLICE_8(:,j)=SWCfull(:,8);
            SWCfull_SLICE_9(:,j)=SWCfull(:,9);
            SWCfull_SLICE_10(:,j)=SWCfull(:,10);
            SWCfull_SLICE_11(:,j)=SWCfull(:,11);
            SWCfull_SLICE_12(:,j)=SWCfull(:,12);
            GDcount_SLICE(j)=GDcount;
            %disp(['Slave ' num2str(j) ' : ' '"Reporting"']);
            %         end
            %[output,CYtic,Runtic,Itic,BMtic]=read_output_yearly(i,j,time,AQflag);
            %^old one, here it reads all the output vars
            %[CY,R]=write_CY_R(output,i,j);
            %             CYstore{time}(i,j)=CY;
            %             Run{time}(i,j)=Run;
            %             I{time}(i,j)=I;
            %             BM{time}(i,j)=BM;


        else
           %disp(['Slave ' num2str(j) ' : ' '"Reporting Else arrival"']);
            %             map(i,j)=0;       %\
            %             CY{time}(i,j)=0;  % | ************************************************
            %             Run{time}(i,j)=0; % | moet herbekeken worden voor als LU erbij komt!!!
            %             I{time}(i,j)=0;   % | ************************************************
            %             BM{time}(i,j)=BM; %/
            %map(jj,ii)=0;       %\
%             CYtemp_SLICE(j)=0;  % | ************************************************
%             Runtemp_SLICE(:,j)=0; % | moet herbekeken worden voor als LU erbij komt!!! LU in rekening bij CN! Voor CY 0 op undeg, best nog een mask doen achteraf via LUmap
%             Itemp_SLICE(j)=0;         % komt niet in actie want --> if landuse(jj+yy1-1,ii+xx1-1) == 1 || landuse(jj+yy1-1,ii+xx1-1)
%             Draintemp_SLICE(j)=0;% | ************************************************
%             BMtemp_SLICE(j)=0; %/
%             SWCtemp_SLICE(j)=0;
%             Pouttemp_SLICE(j)=0;
%             InTree_SLICE(j)=0;
%             %disp(['Slave ' num2str(j) ' : ' '"Reporting Else statement 1"']);
%             STmap_SLICE(j)=0; % anders is STmap te klein wanneer er pixels rechtsonder buiten de ctc wrdn genomen
%             extraP4cell2D_SLICE(j)=0;
%             Evapotemp_SLICE(:,j)=0;
%             Transtemp_SLICE(:,j)=0;
%             Irrtemp_SLICE(:,j)=0;
%             %disp(['Slave ' num2str(j) ' : ' '"Reporting Else statement"']);
        end



        %delete project file from C:\Data_Maarten\software\AquaCrop_batch\LIST\
        %     delete('C:\Data_Maarten\software\AquaCropBatch5\LIST\*.PRO');
       %disp(['Slave ' num2str(j) ' : ' '"Finished work!"']);
    end

   %disp(['"iteration ' num2str(i) ' done."']);

    %% remapping output vars
    subs=sub2ind(size(STmap),rowY_select,colX_select);
    %subs3D=sub2ind(size(Runtemp),rowY_select,colX_select,1:365);
    subs3D=repmat(subs',1,size(Runtemp,3));
    dim3ind=[1:size(Runtemp,3)]-1;
    dim3ind2=repmat(dim3ind,numel(subs),1);
    dim3ind2=dim3ind2(:);
    subs3Dfinal=subs3D'+(dim3ind2.*numel(STmap));

    STmap(subs)=STmap_SLICE';
    Runtemp(subs3Dfinal)=Runtemp_SLICE';
    extraP4cell2D(subs)=extraP4cell2D_SLICE';
    CN_temp(subs)=CN_SLICE;

    CYtemp(subs)=CYtemp_SLICE'; % reden dat deze vars 6x6 zijn, is omdat de Linkse top en zijde 0 krijgen omdat 1e cel is op 2,2, en we later nog extra rij en col bij voegen rechts zijde en onder
    Irrtemp(subs3Dfinal)=Irrtemp_SLICE';
    Itemp(subs)=Itemp_SLICE';
    Draintemp(subs)=Draintemp_SLICE';
    BMtemp(subs)=BMtemp_SLICE';
    SWCtemp(subs)=SWCtemp_SLICE';
    Pouttemp(subs)=Pouttemp_SLICE';
    InTree(subs)=InTree_SLICE';
    Evapotemp(subs3Dfinal)=Evapotemp_SLICE';
    Transtemp(subs3Dfinal)=Transtemp_SLICE';

    StressW_Aer_temp(subs3Dfinal)=StressW_Aer_SLICE';
    StressW_Exp_temp(subs3Dfinal)=StressW_Exp_SLICE';
    StressW_Pol_temp(subs3Dfinal)=StressW_Pol_SLICE';
    StressW_Sen_temp(subs3Dfinal)=StressW_Sen_SLICE';
    StressW_Sto_temp(subs3Dfinal)=StressW_Sto_SLICE';
    StressT_Bio_temp(subs3Dfinal)=StressT_Bio_SLICE';
    StressT_PolC_temp(subs3Dfinal)=StressT_PolC_SLICE';
    StressT_PolH_temp(subs3Dfinal)=StressT_PolH_SLICE';
    SWCfull_temp_1(subs3Dfinal)=SWCfull_SLICE_1';
    SWCfull_temp_2(subs3Dfinal)=SWCfull_SLICE_2';
    SWCfull_temp_3(subs3Dfinal)=SWCfull_SLICE_3';
    SWCfull_temp_4(subs3Dfinal)=SWCfull_SLICE_4';
    SWCfull_temp_5(subs3Dfinal)=SWCfull_SLICE_5';
    SWCfull_temp_6(subs3Dfinal)=SWCfull_SLICE_6';
    SWCfull_temp_7(subs3Dfinal)=SWCfull_SLICE_7';
    SWCfull_temp_8(subs3Dfinal)=SWCfull_SLICE_8';
    SWCfull_temp_9(subs3Dfinal)=SWCfull_SLICE_9';
    SWCfull_temp_10(subs3Dfinal)=SWCfull_SLICE_10';
    SWCfull_temp_11(subs3Dfinal)=SWCfull_SLICE_11';
    SWCfull_temp_12(subs3Dfinal)=SWCfull_SLICE_12';

    extraP4cell_days_temp(subs3Dfinal)=extraP4cell_days_SLICE';

    StressVars.StressW_Aer_temp=StressW_Aer_temp;
    StressVars.StressW_Exp_temp=StressW_Exp_temp;
    StressVars.StressW_Pol_temp=StressW_Pol_temp;
    StressVars.StressW_Sen_temp=StressW_Sen_temp;
    StressVars.StressW_Sto_temp=StressW_Sto_temp;
    StressVars.StressT_Bio_temp=StressT_Bio_temp;
    StressVars.StressT_PolC_temp=StressT_PolC_temp;
    StressVars.StressT_PolH_temp=StressT_PolH_temp;
    SWCVars.SWCfull_temp_1=SWCfull_temp_1;
    SWCVars.SWCfull_temp_2=SWCfull_temp_2;
    SWCVars.SWCfull_temp_3=SWCfull_temp_3;
    SWCVars.SWCfull_temp_4=SWCfull_temp_4;
    SWCVars.SWCfull_temp_5=SWCfull_temp_5;
    SWCVars.SWCfull_temp_6=SWCfull_temp_6;
    SWCVars.SWCfull_temp_7=SWCfull_temp_7;
    SWCVars.SWCfull_temp_8=SWCfull_temp_8;
    SWCVars.SWCfull_temp_9=SWCfull_temp_9;
    SWCVars.SWCfull_temp_10=SWCfull_temp_10;
    SWCVars.SWCfull_temp_11=SWCfull_temp_11;
    SWCVars.SWCfull_temp_12=SWCfull_temp_12;

    GDcount_temp(subs)=GDcount_SLICE';

    colX_select_pool(1:it_num_amount)=[]; % delete the finished it_num from pool
    rowY_select_pool(1:it_num_amount)=[]; % delete the finished it_num from pool

    extraPtemp=Runtemp;
end
nx=xx2-xx1+1;
ny=yy2-yy1;
outerY=zeros(1,nx);
outerX=zeros(1,ny)';

%% check Runnoff - Extra P
Run2D=sum(Runtemp,3);
Col4Run=[1:numel(Run2D(1,:))];
Col4Run=repmat(Col4Run,numel(Run2D(:,1)),1);
Row4Run=[1:numel(Run2D(:,1))]';
Row4Run=repmat(Row4Run,1,numel(Run2D(1,:)));

%buitenste cellen geven aan zichzelf (zodat ze niet naar buiten grenzen
%afvoeren, maakt toch niet uit want voor deze cellen werd geen berekening
%gemaakt boven
rowmatrixcheck=rowmatrix;
colmatrixcheck=colmatrix;
rowmatrixcheck([1 end],:)=0;
rowmatrixcheck(:,[1 end])=0;
colmatrixcheck([1 end],:)=0;
colmatrixcheck(:,[1 end])=0;

RuntempArr=Run2D(:);
Col4RunArr=Col4Run(:);
Row4RunArr=Row4Run(:);
rowmatrixArr=rowmatrixcheck(:);
colmatrixArr=colmatrixcheck(:);
ReceiversRow=Row4RunArr+rowmatrixArr;
ReceiversCol=Col4RunArr+colmatrixArr;

Run2receivingIndx=sub2ind(size(Runtemp),ReceiversRow,ReceiversCol);
[uniqueRec,ia,ic]=unique(Run2receivingIndx);
for iii=1:numel(ia)
    SumUnique(iii)=sum(RuntempArr(ic==iii));
end
Run2Receiving=zeros(size(Runtemp));
% Run2Receiving(Run2receivingIndx)=RuntempArr; % make sure
Run2Receiving(uniqueRec)=SumUnique;
Run2Rec2D=sum(Run2Receiving,3);
DiffExtraPRun=Run2Rec2D-extraP4cell2D;
check_rout_P_extraP=DiffExtraPRun(2:end-1,2:end-1);

%%

% CYtemp=[CYtemp,outerX];%Runoff en SWC(niet waar voor SWC?) hier niet bij want die krijgen al extra rij(+1) + kolom(+1) bij hun initializatie ih begin vd functie
% CYtemp=[CYtemp;outerY];
% % Irrtemp=[Irrtemp,outerX];
% % Irrtemp=[Irrtemp;outerY];
% Itemp=[Itemp,outerX];
% Itemp=[Itemp;outerY];
% Draintemp=[Draintemp,outerX];
% Draintemp=[Draintemp;outerY];
% BMtemp=[BMtemp,outerX];
% BMtemp=[BMtemp;outerY];
% STmap=[STmap,outerX];
% STmap=[STmap;outerY];
% Pouttemp=[Pouttemp,outerX];
% Pouttemp=[Pouttemp;outerY];
% InTree=[InTree,outerX];
% InTree=[InTree;outerY];
% extraP4cell2D=[extraP4cell2D,outerX];
% extraP4cell2D=[extraP4cell2D;outerY];
% % WatBalSum=sum(sum(WatBal));
% % SWCtemp=[SWCtemp,outerX];
% % SWCtemp=[SWCtemp;outerY];
end