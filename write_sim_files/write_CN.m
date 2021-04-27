function [ CN ] = write_CN(landuse, ST,lithology,jj,ii,xx1,yy1,datapath)

%% CN values from table
CN_fallow=[77 86 91 94;76 85 90 93;74 83 88 90];
CN_crops=[72 81 88 91;67 78	85	89];
CN_wood=[45	66 77 83;36 60 73 79;30 55 70 77]; %1= <50% ground cover, 2= 50 to 70% GC, 3= > 70% GC
CN_herba=[80 80 87 93;71 71 81 89;62 62 74 85];
CN_juniper=[75 75 85 89;58 58 73 80;41 41 61 71]; % 1= <30% GC, 2= 30-70% GC, 3= >70% GC

%% loading ST_stone and stone_hydrocond relationships
ST_Stone=load([datapath 'ST_stoneDATA.mat']);
Cond_Stone=load([datapath 'cond_stoneDATA.mat']);

%% Stoniness bepalen
lithval=lithology(jj+yy1-1,ii+xx1-1);
if lithval == 1
    stoninessval=interp1(ST_Stone.ytickval,flipud(ST_Stone.Coll_Stenigheid),ST/10,'spline');
elseif lithval == 2
    stoninessval=interp1(ST_Stone.ytickval,flipud(ST_Stone.Kalk_Stenigheid),ST/10,'spline');
elseif lithval == 3
    stoninessval=interp1(ST_Stone.ytickval,flipud(ST_Stone.Oph_Stenigheid),ST/10,'spline');
elseif lithval == 4
    stoninessval=interp1(ST_Stone.ytickval,flipud(ST_Stone.Con_Stenigheid),ST/10,'spline');
end

if stoninessval < 0
    stoninessval = 0;
end
%% HydroCond bepalen
if lithval == 1 || lithval == 2
    HydroCond=interp1(Cond_Stone.Stoniness,Cond_Stone.Loam,stoninessval,'spline');
elseif lithval == 3 || lithval == 4
    HydroCond=interp1(Cond_Stone.Stoniness,Cond_Stone.SandyLoam,stoninessval,'spline');
end

%% HSG bepalen

if ST < 500
    HSG=4; %D
elseif HydroCond > 1 && HydroCond <= 10
    HSG=3; %C
elseif HydroCond > 10 && HydroCond <= 40
    HSG=2; %B
elseif HydroCond > 40
    HSG=1; %A
end



%% obv LU CN bepalen
if landuse(jj+yy1-1,ii+xx1-1) ==1 %degraded
    CN=CN_crops(2,HSG); % 2 is voor 20% of meer residue on the ground
else % undegraded
    CN=CN_juniper(2,HSG);
end







