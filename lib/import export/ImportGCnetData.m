function [data] = ImportGCnetData(station,tilt_correct)

%% Loading GCnet metadata

filename = '..\AWS_Processing\Input\GCnet\Gc-net_documentation_Nov_10_2000.csv';
delimiter = ';';
startRow = 2;
formatSpec = '%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,4,5]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\.]*)+[\,]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\.]*)*[\,]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers=='.');
                thousandsRegExp = '^\d+?(\.\d{3})*\,{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, '.', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = strrep(numbers, '.', '');
                numbers = strrep(numbers, ',', '.');
                numbers = textscan(numbers, '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

rawNumericColumns = raw(:, [1,3,4,5]);
rawCellColumns = raw(:, 2);

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

overview = table;
overview.ID = cell2mat(rawNumericColumns(:, 1));
overview.StationName = rawCellColumns(:, 1);
overview.Northing = cell2mat(rawNumericColumns(:, 2));
overview.Easting = cell2mat(rawNumericColumns(:, 3));
overview.Elevation = cell2mat(rawNumericColumns(:, 4));

clearvars delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

%% finding station ID
switch station
    case 'CP1'
        station_num = find(strcmp('Crawford Pt.',overview.StationName));
    case 'SouthDome'    
        station_num = find(strcmp('South Dome',overview.StationName));
    case 'NEEM'    
        station_num = 23;
    otherwise
        station_num = find(strcmp(station,overview.StationName));
end

filename = dir(sprintf('../AWS_Processing/Input/GCnet/20190501_jaws/%02i*',station_num));
if size(filename,1) == 0
    filename = dir(sprintf('../AWS_Processing/Input/GCnet/10102018_jaws/%02i*',station_num));
%     if size(filename,1) == 0
%         error('No data file found.')
%     elseif size(filename,1) >1
%         error('Several data files found.')
%     end
elseif size(filename,1) >1
    error('Several data files found.')
else
end
 
%% Loading jaws file

if size(filename,1)~=0
    disp('Reading:')
    disp(filename)
    filename=[filename.folder '\' filename.name];
    file_info = ncinfo(filename);
    list_var = {'Year',                 'year',;...
        'JulianTime',                   'julian_decimal_time';...
        'ShortwaveRadiationDownWm2',    'fsds';...
        'ShortwaveRadiationUpWm2',      'fsus';...
        'NetRadiationWm2',              'fsns';...
        'AirTemperature1C',             'ta_tc1';...
        'AirTemperature2C',             'ta_tc2';...
        'AirTemperature3C',             'ta_cs1';...
        'AirTemperature4C',             'ta_cs2';...
        'RelativeHumidity1Perc',        'rh1';...
        'RelativeHumidity2Perc',        'rh2';...
        'WindSpeed1ms',                 'wspd1';...
        'WindSpeed2ms',                 'wspd2';...
        'WindDirection1deg',            'wdir1';...
        'WindDirection2deg',            'wdir2';...
        'AirPressurehPa',               'ps';...
        'SnowHeight1m',                 'snh1';...
        'SnowHeight2m',                 'snh2';...
        'IceTemperature1C',             'tsn1';...
        'IceTemperature2C',             'tsn2';...
        'IceTemperature3C',             'tsn3';...
        'IceTemperature4C',             'tsn4';...
        'IceTemperature5C',             'tsn5';...
        'IceTemperature6C',             'tsn6';...
        'IceTemperature7C',             'tsn7';...
        'IceTemperature8C',             'tsn8';...
        'IceTemperature9C',             'tsn9';...
        'IceTemperature10C',            'tsn10';...
        'Windspeed2mms',                'wpd_2m';...
        'Windspeed10mms',               'wspd_10m';...
        'WindSensorHeight1m',           'wind_sensor_height_1'; ...
        'WindSensorHeight2m',           'wind_sensor_height_2';...
        'Albedo',                       'alb';...
        'ZenithAngledeg',               'zenith_angle';...
        'SZA_jaws',                     'sza'};
    clearvars data

    data = table;
    vars = {file_info.Variables.Name};
    for i = 1:length(vars)
        aux= strcmp(list_var(:,2),vars{i});
        Index = find(aux);
        if ~isempty(Index)
            data.(list_var{Index,1}) = double(ncread(filename,list_var{Index,2}));
%             disp([list_var{Index,1} '  ' list_var{Index,2} '  ' vars{i}])
%             pause
        else
            try data.(vars{i}) = double(ncread(filename,vars{i}));
            catch me
%                 warning('Could not extract %s',vars{i});
            end
        end
    end

    for i = 1:length(list_var(:,2))
        aux= strcmp(vars(:), list_var(i,2));
        Index = find(aux);
        if ~isempty(Index)
        else
            data.(list_var{i,2})(1) = NaN;
            data.(list_var{i,2})(:) = NaN;
            fprintf('%s is missing from file\n',list_var{i,2});
        end
    end

    % Matlab tie stamp
    data.time = datenum(data.Year,0,data.JulianTime+1);
    
    % using tilt corrected values if asked
    if strcmp(tilt_correct,'yes')
        ind_replace = data.fsds_adjusted~=0;
        data.ShortwaveRadiationDownWm2(ind_replace) = data.fsds_adjusted(ind_replace);
        ind_replace = data.fsus_adjusted~=0;
        data.ShortwaveRadiationUpWm2(ind_replace) = data.fsus_adjusted(ind_replace);
    end
    % removing duplicates
    ind_remove = find(data.time(2:end)-data.time(1:end-1)<=0);
    % disp('Duplicate timesteps:')
    % disp(ind_remove)
    data(ind_remove,:)=[];

    time_2 = datevec(data.time);
    data.MonthOfTheYear = time_2(:,2);
    data.DayOfTheMonth = time_2(:,3);
    data.HourOfTheDay = time_2(:,4);
    data.DayOfTheYear = floor(data.JulianTime+1);
    data.JulianTime = floor(data.JulianTime)+ data.HourOfTheDay/24;
    clearvars time_2
    
    % Extra files from CP1
    if strcmp(station,'CP1')
    
        % Loading first annexe data file
        % The data file that was released by K. Steffen did not contain part of
        % the data. The missing period was communicated later in an 'annex'
        % file. Here we upload it and merge it with the main dataset for CP1.
        data_2 = LoadExtraFileCP1();
        data_new = table;
        data_new.time = unique(union(data.time,data_2.time));
        clearvars ind_neg ind_pos ind_error time_diff temp ind_ok ind_uni

        
            for i = 1:size(data,2)
                if and(~strcmp('time',data.Properties.VariableNames{i}),...
                        isempty(strfind(data.Properties.VariableNames{i},'Origin')))
                    ind_1 = ismember(data_new.time, data.time);
                    data_new.(data.Properties.VariableNames{i}) = NaN(size(data_new.time));
                    data_new.(data.Properties.VariableNames{i})(ind_1) = data.(data.Properties.VariableNames{i});
                    if sum(strcmp((data.Properties.VariableNames{i}),data_2.Properties.VariableNames))>0
                        ind_2 = ismember(data_new.time, data_2.time);
                        data_new.(data.Properties.VariableNames{i})(ind_2) = data_2.(data.Properties.VariableNames{i});
                    end
                end
            end

        data = data_new;
        clearvars data_2 data_new ind_1 ind_2 time_2
    end
  
    data.TemperatureSensorHeight1m = data.WindSensorHeight1m;
    data.TemperatureSensorHeight2m = data.WindSensorHeight2m;
    data.HumiditySensorHeight1m = data.WindSensorHeight1m;
    data.HumiditySensorHeight2m = data.WindSensorHeight2m;
    
    % temperature to degC
    for i = 1:length(data.Properties.VariableNames)
        if ~isempty(strfind(data.Properties.VariableNames{i},'emperatur'))
            data.(data.Properties.VariableNames{i}) = ...
               data.(data.Properties.VariableNames{i}) -273.15;
        end
    end
    
    % pressure to hPa
    data.AirPressurehPa = data.AirPressurehPa/100;
else
    warning('Using text AWS data file')
    [data] = ImportGCnetData_old(station);
end

end
