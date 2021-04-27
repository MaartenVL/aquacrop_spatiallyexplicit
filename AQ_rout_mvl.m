function [Runtemp,CYtemp,BMtemp,Irrtemp,Itemp,STmap,Draintemp,SWCtemp,extraP4cellSum,Pouttemp,InTree,extraP4cell2D,Evapotemp,Transtemp,StressVars,SWCfull_temp,GDcount_temp,CN_temp,extraP4cell_days_temp,check_rout_P_extraP,SWCtemp_INPUTout,Imaxtemp] = AQ_rout_mvl(time, xx1, xx2, yy1, yy2,nrow,ncol,dates_mvl,climateOS,landuse,lithology,sedexport,marshmask,rowY,colX,rowmatrix,colmatrix,outerPix,SWCtemp,testmodeSTshort,tree4aq,CO2atTime,ParamStruct,datapath,AQpath,runoffcalc,slopemap,swc_param,Pprev,swc_paramT)
% reading crop file
CropReadFile=load([datapath 'cropfile_cal4paleo_v4_5aug.mat']); % new, calibrated values of winter wheat jackknife v4
CropReadFileTrees=load([datapath 'cropfile_cal4paleo4trees_v4_5aug.mat']); % new, calibrated values of winter wheat jackknife v4
CropReadFile.CropRead=CropReadFile.CropwheatCal4Paleo;

% lith 4 SWC calc
lithSand=lithology;
lithSand(lithSand==1)=0.4;
lithSand(lithSand==2)=0.3;
lithSand(lithSand==3)=0.7;
lithSand(lithSand==4)=0.5;
lithClay=lithology;
lithClay(lithClay==1)=0.2;
lithClay(lithClay==2)=0.2;
lithClay(lithClay==3)=0.09;
lithClay(lithClay==4)=0.15;

%reading mask file

statusMrow=[1,1,1;0,0,0;-1,-1,-1];
statusMcol=[1,0,-1;1,0,-1;1,0,-1];
extraPtemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
extraP4cellSum=0;
Runtemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
Evapotemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
Transtemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
Irrtemp=zeros(yy2-yy1+1,xx2-xx1+1,365);
SWCfull_temp=zeros(yy2-yy1+1,xx2-xx1+1);
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
SWCtemp_INPUTout=zeros(yy2-yy1+1,xx2-xx1+1);
Imaxtemp=zeros(yy2-yy1+1,xx2-xx1+1);

%if time == 1
    SWCtemp=zeros(yy2-yy1+1,xx2-xx1+1);
%end
%% Preparation parfor loop
% [buurrow,buurcol,rowmatrix,colmatrix,runoffcalc] = routing(dtm,landuse,nrow,ncol,dtm2);
runoffcalcmasked=runoffcalc;
runoffcalcmasked([1 end],:)=0;
runoffcalcmasked(:,[1 end])=0;
%runoffcalcmasked(lithology(yy1:yy2,xx1:xx2)<1)=0; % niet nodig! verknoeit
%volgorde van routing, verkeerde outerpix erafgeknipt!
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

    SWCfull_SLICE=zeros(1,it_num_amount);

%     GDcount_SLICE=zeros(1,it_num_amount);
    CN_SLICE=zeros(1,it_num_amount); 
    extraP4cell_days_SLICE=zeros(365,it_num_amount);
    SWCtemp_INPUT_SLICE=zeros(1,it_num_amount);
    Imaxtemp_SLICE=zeros(1,it_num_amount);

  %parpool(100)
    for j=1:it_num_amount
        
         %j
        w.ProcessId=time;
        %w=getCurrentWorker;
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
            [ST]=get_soil_rout(time,ii,jj,xx1,yy1,sedexport,marshmask,testmodeSTshort,datapath); %TO DO: fix time problem.% TIME SOIL =/= TIME RAIN. Soil in 10yearly, climate = yearly % <<--- fixed the last one
            STmap_SLICE(j)=ST;
            CN=write_CN(landuse,ST,lithology,jj,ii,xx1,yy1,datapath);
            Slope_input=get_slope(slopemap,jj,ii,xx1,yy1);
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
            

            [climateOS_INPUT,Pout,climateOS_Pori]=write_rain_rout365_mvl(climateOS,time,extraP4cell);

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
            
            %if  time == 1; % enkel in de 1e run geen inicon file
                %SWCtemp_INPUT=0.15;%voor SCL is dit de start FC=0.44
                %             SWCtemp=write_IniCon_rout_mvl(ST,ii,jj,SWCtemp); %SWCtemp na model run wordt doorgegeven naar volgend jaar, resultaat hiervan lezen we hier
            %else
            %disp('time=/=1');
                %SWCtemp_INPUT=SWCtemp(jj,ii);
                %end
                if time ==1
                    SWCtemp_INPUT=0.15;
                else
                    if landuse(jj+yy1-1,ii+xx1-1) == 1
                        %SWCtemp_INPUT=swc_param(1)+(swc_param(2)*ST*1000)+(swc_param(3)*Pprev)+(swc_param(4)*runoffcalcmasked(jj-1,ii-1))+(swc_param(5)*double(lithology(jj+yy1-1,ii+xx1-1)))+(swc_param(6)*slopemap(jj+yy1-1,ii+xx1-1));
                        SWCtemp_INPUT=swc_param(1)+(swc_param(2)*ST*1000)+(swc_param(3)*Pprev)+(swc_param(4)*runoffcalcmasked(jj-1,ii-1))+(swc_param(5)*double(lithSand(jj+yy1-1,ii+xx1-1)))+(swc_param(6)*double(lithClay(jj+yy1-1,ii+xx1-1)))+(swc_param(7)*slopemap(jj+yy1-1,ii+xx1-1));
                    else
                        SWCtemp_INPUT=swc_paramT(1)+(swc_paramT(2)*ST*1000)+(swc_paramT(3)*Pprev)+(swc_paramT(4)*runoffcalcmasked(jj-1,ii-1))+(swc_paramT(5)*double(lithSand(jj+yy1-1,ii+xx1-1)))+(swc_paramT(6)*double(lithClay(jj+yy1-1,ii+xx1-1)))+(swc_paramT(7)*slopemap(jj+yy1-1,ii+xx1-1));
                    end
                    if SWCtemp_INPUT<=0 % ik denk dat AQ crasht als SWC =0, dus opheffen met 0.01 (wat ook zeer weinig is)
                        SWCtemp_INPUT=0.01;
                    end
                end

            %Run AQ
            [BM,CY,Irriga,Infill,Drain,Runoff,ET,SWC,SWCfull,Evapo,Trans,StressW,StressT,GDcount,ImaxTot]=AquaCropOS_RUN(climateOS_INPUT,dates_mvl,SoilInfo,CO2atTime,SWCtemp_INPUT,AQpath,j,w,CropReadFile,CropReadFileTrees,LUval,Slope_input,climateOS_Pori,extraP4cell);


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
            SWCfull_SLICE(j)=SWCfull;
            SWCtemp_INPUT_SLICE(j)=SWCtemp_INPUT;
            Imaxtemp_SLICE(j)=ImaxTot;

            
        else

        end

    end

    
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
    SWCtemp_INPUTout(subs)=SWCtemp_INPUT_SLICE';
    
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
    Imaxtemp(subs)=Imaxtemp_SLICE';


    SWCfull_temp(subs)=SWCfull_SLICE';
    
    extraP4cell_days_temp(subs3Dfinal)=extraP4cell_days_SLICE';
    

    
StressVars=0; % MEMORY

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


end