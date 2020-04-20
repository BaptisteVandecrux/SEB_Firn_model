function [time, year, day, hour, pres,...
    T1, T2, z_T1, z_T2, o_T1,o_T2, ...
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, ...
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_out, c] = ...
    ExtractAWSData(c)
%ExtractAWSData: Extract the meteorological data from a csv file located in
%the ./Input folder.
%
% Author: Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================

c.dt_obs  = 3600;   % observational time step in seconds
filename = c.InputAWSFile;

%% ================= Extracting weather data ==============================
%  header
delimiter = '\t';
endRow = 1;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow, 'Delimiter', delimiter,...
    'ReturnOnError', false);
fclose(fileID);
for i = 1:length(dataArray)
    temp = dataArray{i};
    header{i} = temp{1};
end
clearvars endRow formatSpec fileID  ans dataArray

count = length(header);
while isempty(header{count})
    header(count) = [];
    count = count -1;
end
    
% data
startRow = 2;
formatSpec = strcat(repmat('%f', 1,length(header)),'%[^\n\r]');
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
    'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
dataArray (:,(length(header)+1):end) = [];

c.rows = length(dataArray{1});

data_out=table;
for i=1:length(header)
    data_out.(header{i}) = dataArray{:,i};
end
    data_out = standardizeMissing(data_out,-999);
 
% time
year = data_out.Year;
hour = data_out.HourOfDayUTC;
day = data_out.DayOfYear;
% leap years    
time = year + (day + hour/24)/365;
leapyear = find(year/4 == floor(year/4));
if sum(leapyear) >0
    time(leapyear) = year(leapyear)+(day(leapyear)+hour(leapyear)/24.)/366;
end    

%% ========= selecting prescribed period out of the datafile =============
% Here, if param.year is defined, we keep only the melt-year we want
% melt season starts 1st apr until 31 mar of next year. If records starts
% later, takes the begining of the record for start date.
if c.year(1) ~= 0
    if or(c.year>time(end),c.year<time(1))
        error('Choose year within weather data period.')
    end
    year_start = c.year(1);
    if length(c.year)==1
        year_end =year_start+1;
    else
        year_end = c.year(end);
    end
    day_start = datenum(sprintf('01-Apr-%s',c.year))-datenum(sprintf('01-Jan-%s',c.year));
    time_start = year_start + (day_start)/365;
    if year_start/4 == floor(year_start/4)
        time_start = year_start + (day_start)/366;
    end
    if time(1)> time_start
        time_start = time(1);
        disp('Could not start 1st of April. Starting at begining of record.')
    end
    time_end = year_end + (day_start)/365;

    ind = and(time >= time_start , time < time_end);
    data_out = data_out(ind,:);
end

%% ======= interpolation of missing values and basic filtering ============
% if missing data is present for in important fields at the begining or the
% end of the record, then these periods are removed.
data_out = FillLastGaps(data_out,c);

%% ========= Renaming variables ==============
c.M = size(data_out,1);
c.rows = size(data_out,1);

[time, year, day, hour, pres,...
    T1, T2, z_T1, z_T2, o_T1,o_T2, ...
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, ...
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs] = ...
    RenamingVariables(data_out, c);

%% Loading alternative thermistor string for FA site
if strcmp(c.station,'Miege')
    [T_ice_obs, depth_thermistor] = ...
        LoadThermFA13(data_out.time, data_out.SurfaceHeightm);
end

end
