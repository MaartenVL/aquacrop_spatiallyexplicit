function [ dates_mvl,dates4climOS ] = get_time_period
Year_start=2000;
Month_start=9;
Day_start=1;
Year_end=2001;
Month_end=8;
Day_end=31;

if Month_start < 10 && Day_start < 10
    startdate=[num2str(Year_start) '-0' num2str(Month_start) '-0' num2str(Day_start)];
elseif Month_start < 10
    startdate=[num2str(Year_start) '-0' num2str(Month_start) '-' num2str(Day_start)];
elseif Day_start < 10
    startdate=[num2str(Year_start) '-' num2str(Month_start) '-0' num2str(Day_start)];
end

if Month_end < 10 && Day_end < 10
    enddate=[num2str(Year_end) '-0' num2str(Month_end) '-0' num2str(Day_end)];
elseif Month_end < 10
    enddate=[num2str(Year_end) '-0' num2str(Month_end) '-' num2str(Day_end)];
elseif Day_end < 10
    enddate=[num2str(Year_end) '-' num2str(Month_end) '-0' num2str(Day_end)];
end

dates_mvl={startdate;enddate};
dates_mvl_numeric=[Year_start,Month_start,Day_start;Year_end,Month_end,Day_end]

% days 4 climateOS file
month=[31 28 31 30 31 30 31 31 30 31 30 31];
monthcum=cumsum(month);
monthcum=[0 monthcum];
monthcum(end)=[];
for i=1:numel(month)
    DayAr{i}=linspace(1,month(i),month(i));
    MoAr1=zeros(size(DayAr{i}))+i;
    MoAr{i}=MoAr1;
    clear MoAr1
end

Day2years=[DayAr{:} DayAr{:}];
Month2years=[MoAr{:} MoAr{:}];
Year2years=[repmat(2000,1,365) repmat(2001,1,365)];
ind4start=monthcum(Month_start)+Day_start;
ind4end=ind4start+364;
Days=Day2years(ind4start:ind4end);
Months=Month2years(ind4start:ind4end);
Years=Year2years(ind4start:ind4end);

dates4climOS=[Days',Months',Years'];
