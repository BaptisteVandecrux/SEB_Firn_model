function [data_RCM] = LoadRCMData(station,RCM)

% Loading AWS list
opts = delimitedTextImportOptions("NumVariables", 5);
opts.DataLines = [2, Inf];
opts.Delimiter = ";";
opts.VariableNames = ["ID", "Stationname", "Northing", "Easting", "Elevationmasl"];
opts.VariableTypes = ["double", "string", "double", "double", "double"];
opts = setvaropts(opts, 2, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
AWS = readtable("../AWS_Processing/Input/Secondary data/AWS.csv", opts);
clear opts

ind_station = find(strcmp(AWS.Stationname,station));

% loading RCM data at the station
filename = sprintf('../AWS_Processing/Input/Secondary data/%s_3h_AWS_sites.nc',RCM);
if ~exist(filename,'file')
   filename = sprintf('../AWS_Processing/Input/Secondary data/%s_hourly.nc',RCM);
end
finfo = ncinfo(filename);
varname = {finfo.Variables.Name }';
data_RCM = table();
try data_RCM.time = ncread(filename,'time')+datenum(1900,1,1);
catch me 
    data_RCM.time = ncread(filename,'Time')+datenum(1900,1,1);
end

for i = 3: length(varname)
    tmp = ncread(filename,varname{i});
    if size(tmp,2)>1
        tmp = tmp(:,ind_station);
    end

    data_RCM.(varname{i}) = tmp;
end

% Resampling at hourly value
data_RCM = standardizeMissing(data_RCM,-999);
if strcmp(RCM,'RACMO')
    data_RCM(:,11:12) = [];
end
data_RCM_org = data_RCM;
data_RCM = ResampleTable(data_RCM_org,'spline');   
data_RCM2 = ResampleTable(data_RCM_org,'previous');   

var_pr = {'prsn','pr','rf','sf'};
for i = 1:length(var_pr)
    if contains(data_RCM.Properties.VariableNames, var_pr{i})
        data_RCM.(var_pr{i}) = data_RCM.(var_pr{i});
    end
end

% harmonizing things
if ismember('tas',data_RCM.Properties.VariableNames)
    data_RCM.ta2m = data_RCM.tas;
    data_RCM.tas = [];
end

if ismember('tas2m',data_RCM.Properties.VariableNames)
    data_RCM.ta2m = data_RCM.tas2m;
    data_RCM.tas2m = [];
end

if ismember('prsn',data_RCM.Properties.VariableNames)
    data_RCM.sf = data_RCM.prsn;
    data_RCM.rf = data_RCM.pr-data_RCM.prsn;
    data_RCM.prsn = [];
else
    if ismember('pr',data_RCM.Properties.VariableNames)
        data_RCM.sf = data_RCM.pr;
    end
end

if ismember('relhum2m',data_RCM.Properties.VariableNames)
    data_RCM.relhum = data_RCM.relhum2m;
    data_RCM.relhum2m = [];
end

if ismember('rhs',data_RCM.Properties.VariableNames)
    data_RCM.relhum = data_RCM.rhs;
    data_RCM.rhs = [];
end

if ~ismember('ps',data_RCM.Properties.VariableNames)
    if ismember('psl',data_RCM.Properties.VariableNames)
        data_RCM.ps = data_RCM.psl/100;
        data_RCM.psl = [];
    end
end

if ~ismember('relhum',data_RCM.Properties.VariableNames)
    c.T_0 = 273.15;
    c.T_100 = 373.15;
    c.es = 0.622;
    c.es_0 = 6.1071;
    c.es_100 = 1013.246;

    data_RCM.relhum = spechum2relhum(data_RCM.ta2m+c.T_0,...
        data_RCM.ps,data_RCM.q2m*0.001,c);
end

%% Renaming variables
    IndexC = strfind(data_RCM.Properties.VariableNames, 'ta2m');
    ind_var = find(not(cellfun('isempty', IndexC)));
    data_RCM.Properties.VariableNames{ind_var} = 'AirTemperature2C';
    data_RCM.AirTemperature1C = data_RCM.AirTemperature2C;
       
    IndexC = strfind(data_RCM.Properties.VariableNames, 'dswrad');
    if isempty( [IndexC{:}])
            IndexC = strfind(data_RCM.Properties.VariableNames, 'rsds');
    end
    ind_var = find(not(cellfun('isempty', IndexC)));
    data_RCM.Properties.VariableNames{ind_var} = 'ShortwaveRadiationDownWm2';

    IndexC = strfind(data_RCM.Properties.VariableNames, 'ps');
    ind_var = find(not(cellfun('isempty', IndexC)));
    data_RCM.Properties.VariableNames{ind_var} = 'AirPressurehPa';

    IndexC = strfind(data_RCM.Properties.VariableNames, 'relhum');
    ind_var = find(not(cellfun('isempty', IndexC)));
    data_RCM.Properties.VariableNames{ind_var} = 'RelativeHumidity2Perc';
    data_RCM.RelativeHumidity1Perc = data_RCM.RelativeHumidity2Perc;

    IndexC = strfind(data_RCM.Properties.VariableNames, 'wind10');
    if isempty( [IndexC{:}])
            IndexC = strfind(data_RCM.Properties.VariableNames, 'sfcWind');
    end
    ind_var = find(not(cellfun('isempty', IndexC)));
    data_RCM.Properties.VariableNames{ind_var} = 'WindSpeed2ms';
    data_RCM.WindSpeed1ms = data_RCM.WindSpeed2ms;
    
    IndexC = strfind(data_RCM.Properties.VariableNames, 'dlwrad');
    if isempty( [IndexC{:}])
            IndexC = strfind(data_RCM.Properties.VariableNames, 'rlds');
    end
    ind_var = find(not(cellfun('isempty', IndexC)));
    data_RCM.Properties.VariableNames{ind_var} = 'LongwaveRadiationDownWm2';
    
    % assigning heights
    data_RCM.WindSensorHeight1m = 10 * ones(size(data_RCM,1),1);
    data_RCM.WindSensorHeight2m = 10 * ones(size(data_RCM,1),1);
    data_RCM.TemperatureSensorHeight1m = 2 * ones(size(data_RCM,1),1);
    data_RCM.TemperatureSensorHeight2m  = 2 * ones(size(data_RCM,1),1);
    data_RCM.HumiditySensorHeight1m = 2 * ones(size(data_RCM,1),1);
    data_RCM.HumiditySensorHeight2m = 2 * ones(size(data_RCM,1),1);

    data_RCM.Snowfallmweq = data_RCM.sf;
    if ismember('rf',data_RCM.Properties.VariableNames)      
        data_RCM.Rainfallmweq = data_RCM.rf;
        data_RCM.rf = [];
    else
        data_RCM.Rainfallmweq = data_RCM.Snowfallmweq*0;
    end
    
    data_RCM.Snowfallmweq = data_RCM.Snowfallmweq/1000;
    data_RCM.Rainfallmweq = data_RCM.Rainfallmweq/1000;

    if ismember('Rainfallmweq',data_RCM.Properties.VariableNames)
        figure
        plot(data_RCM.time,cumsum(data_RCM.Rainfallmweq))
        title(sprintf('Cumulated rain at %s',station))
        datetick('x')
        axis tight
        ylabel('m weq')   
    end
end