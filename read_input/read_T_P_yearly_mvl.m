function [Tmin,Tmax,T,P,Pprev] = read_T_P_yearly_mvl(time,dtime,T_HR,P_HR_rand,Tmax_in,Tmin_in)

%% aanmaken T en P slices for each year
T=T_HR(time*10,:); %meanT for 10 years
T=circshift(T,[1 4]); %need to have September be the first month from which we read Temp since this is also first month of simulation
Tmax=Tmax_in(time*10,:);
Tmin=Tmin_in(time*10,:);
Tmax= circshift(Tmax,[1 4]); %need to have September be the first month from which we read Temp since this is also first month of simulation
Tmin= circshift(Tmin,[1 4]); %need to have September be the first month from which we read Temp since this is also first month of simulation
%P=mean(P_HR_rand((time_actual*10)-9:time_actual*10,:)); %mean P for 10 years niet doen,
%kies gewoon 1 jaar, anders gaan je neerslag 0 dagen verloren...
P=P_HR_rand(time*10,:);
P=circshift(P,[1 30+30+31+31]); %need to have September be the first month from which we read Precip since this is also first month of simulation

if time~=1
    Pprev=sum(P_HR_rand((time*10)-1,:),2); %year before
    %Pprev=sum(sum(P_HR_rand(((time-dtime)*10):((time*10)-1),:),2))/20; % average past 20 yr...via dtime
else
    Pprev=0;
end

%% testing Pprev - T tested manually...
% time=[1:2:400];
% P1=sum(P_HR_rand(time.*10,:),2);
% Pprev=sum(P_HR_rand((time.*10)-1,:),2);
% Psum=sum(P_HR_rand,2); %look P1 and Pprec up in this one, it matches...
end

