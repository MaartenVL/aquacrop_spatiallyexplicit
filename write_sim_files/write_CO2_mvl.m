function [ CO2atTime ] = write_CO2_mvl(time,dtime,datapath)
%% CO2 import
filename = [datapath 'CO2 NOAA.txt'];
delimiter = '\t';
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
CO2 = dataArray{:, 1};
year = dataArray{:, 2};
clearvars filename delimiter formatSpec fileID dataArray ans;

% year=year+1; %zodat 0bp model time 1 wordt
% year=flipud(year);
% diff=abs(year-time);
% [val,ind]=min(diff);
% CO2atTime=CO2(ind);

% Convert time 2 BP
%time_full=fliplr([1:dtime*10:4000]);
%BPyear=time_full(time);
BPyear=4000-(time*10);
diff=abs(year-BPyear);
[val,ind]=min(diff);
indarray=[ind+1 ind ind-1];
if indarray(1) > numel(CO2)
    indarray(1)=indarray(2)
end
if indarray(3) < 1
    indarray(3)=indarray(2);
end
CO2array=CO2(indarray);
CO2atTime=mean(CO2array);

%% testing code
% clear all
% datapath='/data/leuven/318/vsc31850/AQ_on_HPC_parallel_v3/AQ_OSmvl2/input_mvl_files/';
% filename = [datapath 'CO2 NOAA.txt'];
% delimiter = '\t';
% formatSpec = '%f%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% CO2 = dataArray{:, 1};
% year = dataArray{:, 2};
% clearvars filename delimiter formatSpec fileID dataArray ans;
% 
% %  dtime=20;
% time=[1:2:400];
% % time_full=fliplr([1:dtime*10:4000]);
% % BPyear=time_full(time);
% 
% BPyear=4000-(time.*10);
% for i=1:numel(BPyear)
% diff=abs(year-BPyear(i));
% [val(i),ind(i)]=min(diff);
% CO2atTime(i)=CO2(ind(i));
% 
% indarray=[ind(i)+1 ind(i) ind(i)-1];
% if indarray(1) > numel(CO2)
%     indarray(1)=indarray(2)
% end
% if indarray(3) < 1
%     indarray(3)=indarray(2);
% end
% CO2array=CO2(indarray);
% CO2atTime_mean(i)=mean(CO2array);
% end
% figure
% plot(CO2atTime_mean);
% hold on
% plot(CO2atTime,'r');
end
