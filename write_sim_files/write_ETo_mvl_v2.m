function [EToday,Rn,Tday,TmaxDay,TminDay] = write_ETo_mvl_v2(Tmin,Tmax,T)
%%ETo calc
%[STEP 1] Cummulative days in the year
%-------------------------------------
% T=circshift(T,[1 -4]);
% Tmin=circshift(Tmin,[1 -4]);
% Tmax=circshift(Tmax,[1 -4]);



month=[31 28 31 30 31 30 31 31 30 31 30 31];
month2=cumsum(month);

timevect=linspace(1,365,365);
month3=month2-(month./2); % middle of each month, cummulative days, xvector of linear interpolation
Tday=interp1(month3,T, timevect,'spline');
TmaxDay=interp1(month3,Tmax, timevect,'spline');
TminDay=interp1(month3,Tmin, timevect,'spline');

%[STEP 2] Creating empty variables
%----------------------------------
d=zeros(365);
ws=zeros(365);
Ra=zeros(365);
EThar=zeros(365);
Dr=zeros(365);
RaAvg=zeros(12);
Rsum=0;

%[STEP 3] Calculating Ra
%-----------------------
i=linspace(1,365,365);
d=0.409.*(sin((2*3.14.*(i)/365)-1.39));
Dr=1+(0.033.*cos(2*3.14.*(i)/365));
lat=0.65;
ws=acos(-tan(lat).*tan(d));
Ra=((24*60)/3.14).*0.0820.*Dr.*((ws.*sin(lat).*sin(d))+(cos(lat).*cos(d).*sin(ws)));

%shifting days to sep start
dayshift=30+30+31+31;
circshift(Ra,[1 dayshift]); %from jan start, to a sep start.

% Using the followong Hargreaves equation, ETo is calculated
%Gavilan et al 2006
ETo=zeros(1,12);

%Allen 1998
G=0;
sigma=4.903*(10^(-9));
delta=4098*(0.6108.*exp((17.27.*Tday)./(Tday+237.3)))./((Tday+237.3).^2);
Tdew=TminDay-2;
ea=0.6108.*exp((17.27.*Tdew)./(Tdew+237.3));
esTmaxDay=0.6108.*exp((17.27.*TmaxDay)./(TmaxDay+237.3));
esTminDay=0.6108.*exp((17.27.*TminDay)./(TminDay+237.3));
es=(esTmaxDay+esTminDay)/2;
Rs=0.16.*sqrt(TmaxDay-TminDay).*Ra;
Rns=(1-0.23).*Rs;
Rso=(0.75+((2*(10^(-5))*1200)))*Ra;
Rnl=sigma*((((TmaxDay+273.16).^4)+((TminDay+273.16).^4))/2).*(0.34-(0.14.*sqrt(ea))).*((1.35*(Rs./Rso))-0.35);
Rn=Rns-Rnl;
u2=2;
epsi=1.013*(10^(-3))*101.3/0.622*2.45;
ETo=(0.408.*delta.*(Rn-G)+(epsi*(900./(Tday+273)).*u2.*(es-ea)))./(delta+(epsi*(1+(0.34*u2))));

EToday=ETo;

end


