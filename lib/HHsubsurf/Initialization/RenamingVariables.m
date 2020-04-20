function [time, year, day, hour, pres,...
    T1, T2, z_T1, z_T2, o_T1,o_T2, ...
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, ...
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,...
    SRin,SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs] = ...
    RenamingVariables(data_out,c)

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
    
    % temperature, humidity and wind speed
    if sum(strcmp(data_out.Properties.VariableNames,'AirTemperatureC'))
        disp('Only one level was detected on the weather station.')
        T2 = data_out.AirTemperatureC;
        T1 = NaN(size(T2));
        RH2 = data_out.RelativeHumidity;
        RH1 = NaN(size(T2));
        WS2 = data_out.WindSpeedms;
        WS1 = NaN(size(T1));

        o_T1 = NaN(size(T1));
        o_RH1 = NaN(size(T1));
        o_WS1 = NaN(size(T1));
        z_T1 = NaN(size(T1));
        z_RH1 = NaN(size(T1));
        z_WS1 = NaN(size(T1));

        % assigning height
        if sum(strcmp('HeightWindSpeedm',data_out.Properties.VariableNames))
            z_WS2 = data_out.HeightWindSpeedm;
            z_T2 = data_out.HeightTemperaturem;
            z_RH2 = data_out.HeightHumiditym;
        else
            disp('Assigning measurement height from HeightSensorBoomm_raw field.')
            z_WS2 = data_out.HeightSensorBoomm_raw + 0.4;
            z_T2 = data_out.HeightSensorBoomm_raw - 0.12;
            z_RH2 = data_out.HeightSensorBoomm_raw - 0.12;
        end

        % assigning origin
        if sum(strcmp('WindSpeed1ms_Origin',data_out.Properties.VariableNames))
            o_WS2 = data_out.WindSpeed1ms_Origin;
            o_T2 = data_out.AirTemperature1C_Origin;
            o_RH2 = data_out.RelativeHumidity1_Origin;
        else
            disp('No WindSpeed1ms_Origin field specified')
            o_T2 = zeros(size(T1));
            o_RH2 = zeros(size(T1));
            o_WS2 = zeros(size(T1));
        end
        
    elseif sum(strcmp(data_out.Properties.VariableNames,'AirTemperature1C'))
        disp('Two levels detected on the weather station')
        T1 = data_out.AirTemperature1C;
        T2 = data_out.AirTemperature2C;
        RH1 = data_out.RelativeHumidity1;
        RH2 = data_out.RelativeHumidity2;
        WS1 = data_out.WindSpeed1ms;
        WS2 = data_out.WindSpeed2ms;
        
        o_T1 = data_out.AirTemperature1C_Origin;
        o_T2 = data_out.AirTemperature2C_Origin;
        o_RH1 = data_out.RelativeHumidity1_Origin;
        o_RH2 = data_out.RelativeHumidity2_Origin;
        o_WS1 = data_out.WindSpeed1ms_Origin;
        o_WS2 = data_out.WindSpeed2ms_Origin;

        z_T1 = data_out.HeightTemperature1m;
        z_T2 = data_out.HeightTemperature2m;
        z_RH1 = data_out.HeightHumidity1m;
        z_RH2 = data_out.HeightHumidity2m;
        z_WS1 = data_out.HeightWindSpeed1m;
        z_WS2 = data_out.HeightWindSpeed2m;
    else
        error('Cannot recognize temperature field in weather data file.')
    end
    T1 = T1 + c.T_0;       % Temperature in degrees Kelvin
    T2 = T2 + c.T_0;       % Temperature in degrees Kelvin

    % radiation
    LRin = data_out.LongwaveRadiationDownWm2;
    LRout = data_out.LongwaveRadiationUpWm2;

    if sum(strcmp(data_out.Properties.VariableNames,'ShortwaveRadiationDown_CorWm2'))
        disp('Using ShortwaveRadiation***_CorWm2 field.')
        SRin = data_out.ShortwaveRadiationDown_CorWm2;
        SRout = data_out.ShortwaveRadiationUp_CorWm2;
    else
        SRin = data_out.ShortwaveRadiationDownWm2;
        SRout = data_out.ShortwaveRadiationUpWm2;
    end

    % other variables
    pres = data_out.AirPressurehPa;
    Surface_Height = data_out.SurfaceHeightm;   
    Tsurf_obs = min(c.T_0, ((LRout-(1-c.em)*LRin)/c.em/c.sigma).^0.25);
    Tsurf_obs(or(isnan(LRout),isnan(LRin))) = NaN;
    
    ind = strfind(data_out.Properties.VariableNames,'IceTemperature');
    ind = find(~cellfun('isempty', ind));
    ind2 = strfind(data_out.Properties.VariableNames,'DepthThermistor');
    ind2 = find(~cellfun('isempty', ind2));
    num_therm = length(ind);
    T_ice_obs = NaN(length(T1),num_therm);
    depth_thermistor = NaN(length(T1),num_therm);

    if ~isempty(ind2)
        for i = 1:length(ind)
            T_ice_obs(:,i) = data_out.(data_out.Properties.VariableNames{ind(i)});
            depth_thermistor(:,i) = data_out.(data_out.Properties.VariableNames{ind2(i)});
        end
    end
        

    ind =  (LRout>316);
    if sum(ind)>0
        if c.verbose == 1
        fprintf('Warning: Observed surface temperature higher than 0degC\n')
        end
    %     before = LRout(ind);
    % %     LRout(ind) = LRout(ind) - (20/15 * (T(ind)-c.T_0));
    %     figure
    %     scatter(T(ind),before,'xr')
    % %     hold on
    % %     scatter(T(ind),LRout(ind),'ob')
    % %     legend('before correction', 'after correction')
    %     xlabel('Temperature (deg C)')
    %     ylabel('Outgoing long-wave radiation (W/m^2)')
    %     title('Observed LRout > black body at 0degC')
    end
end