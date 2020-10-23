function [data_AWS,Tsurf_obs, pres, T1, T2, z_T1, z_T2, ...
    o_T1, o_T2,RH1, RH2, z_RH1, z_RH2, ...
    o_RH1, o_RH2,WS1, WS2, z_WS1, z_WS2, ...
    o_WS1, o_WS2, SRin, SRout, LRin, LRout, c] = PrepareForRetMIP(c)

% for the RetMIP forcing, we need to by-pass the weather data forcing
    % and use the prescribed melt and skin temperature directly
    filename = sprintf('./RetMIP/Input files/surface/RetMIP_%s.csv',c.station);
    delimiter = ';';
    startRow = 2;
    formatSpec = '%{dd-MMM-yyyy HH:mm:ss}D%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    data_AWS = table(dataArray{1:end-1}, 'VariableNames', {'time','melt_mmweq','acc_subl_mmweq','Tsurf_K'});
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    % time
    year = data_AWS.time.Year;
    hour = data_AWS.time.Hour;
    DV = datevec(data_AWS.time);
    day = 	datenum(DV(:,1),DV(:,2),DV(:,3))-datenum(DV(:,1),0,0);
    % leap years    
    time = year + (day + hour/24)/365;
    leapyear = find(year/4 == floor(year/4));
    if sum(leapyear) >0
        time(leapyear) = year(leapyear)+(day(leapyear)+hour(leapyear)/24.)/366;
    end
    
    Tsurf_obs = data_AWS.Tsurf_K;
    c.solve_T_surf = 0;
    
    c.dt_obs  = (hour(2) - hour(1)) *3600;   % observational time step in seconds
    c.zdtime = c.dt_obs;
    c.delta_time = c.dt_obs;
    c.M = size(data_AWS,1);
    c.rows = size(data_AWS,1);
    c.rho_snow = ones(c.rows,1)*315;
    
    data_AWS.Snowfallmweq = data_AWS.acc_subl_mmweq/1000;
    
    pres= NaN(c.rows,1);
    T1= NaN(c.rows,1);
    T2= NaN(c.rows,1);
    z_T1= NaN(c.rows,1);
    z_T2= NaN(c.rows,1);
    o_T1= NaN(c.rows,1);
    o_T2= NaN(c.rows,1);
    RH1= NaN(c.rows,1);
    RH2= NaN(c.rows,1);
    z_RH1= NaN(c.rows,1);
    z_RH2= NaN(c.rows,1);
    o_RH1= NaN(c.rows,1);
    o_RH2= NaN(c.rows,1);
    WS1= NaN(c.rows,1);
    WS2= NaN(c.rows,1);
    z_WS1= NaN(c.rows,1);
    z_WS2= NaN(c.rows,1);
    o_WS1= NaN(c.rows,1);
    o_WS2= NaN(c.rows,1);
    SRin= NaN(c.rows,1);
    SRout= NaN(c.rows,1);
    LRin= NaN(c.rows,1);
    LRout= NaN(c.rows,1);
end
