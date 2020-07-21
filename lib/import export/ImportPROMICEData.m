function [data] = ImportPROMICEData(station)
    % Hourly values from the station
        filename = sprintf('../AWS_Processing/Input/PROMICE/%s_hour_v03.txt',station);

        delimiter = ' ';
        startRow = 2;
        formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        fclose(fileID);

        data = table(dataArray{1:end-1}, 'VariableNames', ...
            {'Year','MonthOfYear','DayOfMonth','HourOfDayUTC','DayOfYear',...
            'DayOfCentury','AirPressurehPa','AirTemperature1C','AirTemperature3C',...
            'RelativeHumidity1Perc','SpecificHumiditygkg','WindSpeed1ms','WindDirection1deg',...
            'SensibleHeatFluxWm2','LatentHeatFluxWm2','ShortwaveRadiationDown_unCorWm2',...
            'ShortwaveRadiationDownWm2','ShortwaveRadiationUp_unCorWm2',...
            'ShortwaveRadiationUpWm2','Albedo_theta70d','LongwaveRadiationDownWm2',...
            'LongwaveRadiationUpWm2','CloudCover','SurfaceTemperatureC',...
            'HeightSensorBoomm','HeightStakesm','DepthPressureTransducerm',...
            'DepthPressureTransducer_Corm','IceTemperature1C','IceTemperature2C',...
            'IceTemperature3C','IceTemperature4C','IceTemperature5C',...
            'IceTemperature6C','IceTemperature7C','IceTemperature8C',...
            'TiltToEastd','TiltToNorthd','TimeGPShhmmssUTC','LatitudeGPSdegN',...
            'LongitudeGPSdegW','ElevationGPSm','HorDilOfPrecGPS',...
            'LoggerTemperatureC','FanCurrentmA','BatteryVoltageV'});

        clearvars filename delimiter startRow formatSpec fileID dataArray ans;

        data = standardizeMissing(data,-999);
        data.time = datenum(data.Year,data.MonthOfYear,data.DayOfMonth,data.HourOfDayUTC,0,0);

      
    switch station
        case 'KAN_U'
        %% data specific to KANU:
% during a period, the data from KANU was transmitted but not saved on
% we therefor need to incorporate these daily transmission in the
% timeseries. This was done by Babis so we just take his values as
% secondary values in the gap-filling procedure.
       
%         data_daily = import_daily(sprintf('%s_day_v02.txt',station));
%         data_daily = standardizeMissing(data_daily,-999);
%         data_daily.time = datenum(data_daily.Year,data_daily.MonthOfYear,data_daily.DayOfMonth,0,0,0);
%         
% %         ind_remove = and(data.time>=datenum('08-Sep-2010 23:00:00') , ...
% %             data.time<=datenum('22-Oct-2011 23:00:00'));
% %         data.HeightSensorBoomm(ind_remove) = NaN;
%         
%         % Using values from daily file to replace missing hourly values
%         ind_fromdaily = [];
%         for i = 1:size(data_daily,1)
%             if isnan(data.IceTemperature3C(data.time==data_daily.time(i)))
%                 ind_fromdaily = [ind_fromdaily i];
%                 ind = (data.time==data_daily.time(i));
%                 data.HeightStakesm(ind) ...
%                     = data_daily.HeightStakesm(i);
%                 data.HeightSensorBoomm(ind) ...
%                     = data_daily.HeightSensorBoomm(i);
%                 data.IceTemperature1C(ind) ...
%                     = data_daily.IceTemperature1C(i);
%                 data.IceTemperature2C(ind) ...
%                     = data_daily.IceTemperature2C(i);
%                 data.IceTemperature3C(ind) ...
%                     = data_daily.IceTemperature3C(i);
%                 data.IceTemperature4C(ind) ...
%                     = data_daily.IceTemperature4C(i);
%                 data.IceTemperature5C(ind) ...
%                     = data_daily.IceTemperature5C(i);
%                 data.IceTemperature6C(ind) ...
%                     = data_daily.IceTemperature6C(i);
%                 data.IceTemperature7C(ind) ...
%                     = data_daily.IceTemperature7C(i);
%                 data.IceTemperature8C(ind) ...
%                     = data_daily.IceTemperature8C(i);
%             end
%         end
% 
%         f = figure;
%         plot(data.time,data.HeightSensorBoomm,'LineWidth',2)
%         hold on
%         scatter(data_daily.time(ind_fromdaily), data_daily.HeightSensorBoomm(ind_fromdaily))
%         datetick('x','dd-mm-yyyy', 'keeplimits', 'keepticks')
%         xlabel('Date')
%         ylabel('Height of sensor boom (m)')
%         legend('Hourly data','Daily data used to fill the gaps')
        
        case 'NUK_K'
            RH_wrt_i = data.RelativeHumidity1Perc;
            T = data.AirTemperature3C + 273.15;
            pres = data.AirPressurehPa*100;
            data.RelativeHumidity1Perc = RHice2water(RH_wrt_i,T,pres);
    end

    data.WindSensorHeight1m  = data.HeightSensorBoomm +0.4;
    data.TemperatureSensorHeight1m  = data.HeightSensorBoomm -0.1;
    data.HumiditySensorHeight1m  = data.HeightSensorBoomm -0.1;
    
    varlist1 = {'AirTemperature1C','RelativeHumidity1Perc','WindSpeed1ms','WindDirection1deg',...
        'WindSensorHeight1m','TemperatureSensorHeight1m','HumiditySensorHeight1m'};
    varlist2 = {'AirTemperature2C','RelativeHumidity2Perc','WindSpeed2ms','WindDirection2deg',...
        'WindSensorHeight2m','TemperatureSensorHeight2m','HumiditySensorHeight2m'};
    for k = 1:length(varlist1)
        data.(varlist2{k}) = data.(varlist1{k});  
    end
    data.AirTemperature4C = NaN(size(data.AirTemperature1C));
    data.NetRadiationWm2 = NaN(size(data.RelativeHumidity1Perc));
    if sum(~isnan(data.HeightStakesm))>10
        data.SnowHeight2m = data.HeightStakesm(find(~isnan(data.HeightStakesm),1,'first')) - data.HeightStakesm;
    else
        data.SnowHeight2m = NaN(size(data.RelativeHumidity1Perc));
    end
    data.SnowHeight1m = data.WindSensorHeight1m(find(~isnan(data.WindSensorHeight1m),1,'first')) - data.WindSensorHeight1m;
end
    