function data = TreatAndFilterData(data, opt, station,OutputFolder,vis)
% data processing and filtering

% in the text files missing data is either 999 or -999
% we change them into Matlab's 'NaN'
data = standardizeMissing(data,{999,-999});
% data = ResampleTable(data);

%% Site specific modification
    switch station
        case 'CP1'
            [data] = SpecialTreatmentCP1(data);
            % When generating the CP1 data, the surface height from CP2 is being used.
            % This is the only station where surface height is reconstructed from
            % another location. For all the other station the output from HIRHHAM5 will
            % be used instead
            disp('Recover Surface Height From CP2')
            tic
            data = RecoverSurfaceHeightFromCP2(data,'CP2',OutputFolder,vis);
            toc
        case 'DYE-2'
            DV  = datevec(data.time);  % [N x 6] array

            ind2 = and(DV(:,4)>=21,DV(:,4)<=23);
            ind = and(ind2,data.time>= datenum('01-Jan-2014'));
            data.ShortwaveRadiationDownWm2(ind)=NaN;  
            
            ind1= find(data.time==datenum('15-May-2014'));
            ind2= find(data.time==datenum('22-Apr-2016'));
            data.AirPressurehPa(ind1:ind2) = data.AirPressurehPa(ind1:ind2)-198.8;
                        
            ind1= find(data.time==datenum('20-Oct-2016'));
            data.AirPressurehPa(ind1:end) = data.AirPressurehPa(ind1:end)-198.8;
            data.AirPressurehPa(data.AirPressurehPa>825) = NaN;
        case 'NASA-SE'
            ind1= find(data.time==datenum('01-May-2017'));
            data.AirPressurehPa(ind1:end) = data.AirPressurehPa(ind1:end)+404;
        
        case 'Summit'
          ind1= find(data.time==datenum( '02-May-2017'));
            data.AirPressurehPa(ind1:end) = data.AirPressurehPa(ind1:end)-50;
        case 'NUK_K'
            %measured snow thickness at installation
            data.SnowHeight2m=data.SnowHeight2m+1.56;
            [data] = SpecialTreatmentNUK_K(data);

        case 'NASA-U'
            [data] = SpecialTreatmentNASAU(data);
        case 'Saddle'
            data(find(abs(data.time-datenum('23-May-2013 13:59:57'))<1/24),:) =...
                [];
            ind1= find(data.time==datenum('01-May-2017'));
            data.AirPressurehPa(ind1:end) = data.AirPressurehPa(ind1:end)+404;
        case 'GITS'
            data.AirTemperature3C(...
                data.AirTemperature1C>nanmax(data.AirTemperature2C)*1.05)=NaN;
    end
%% Correcting radiation
% Suggested at some point by Robert
%     disp('correcting Radiation'
% tic
%     data = CorrectingRadiation(data,station,vis);
% toc

%% converting humidity
% The instruments on PROMICE stations or HIRHAM model normally give 
% humidity with regard to water, but GCnet instruments give it with regards
% to ice. So we convert it back to humidity with regard to water for
% compatibility.
if strcmp(opt,'ConvertHumidity')
    RH_wrt_i = data.RelativeHumidity1Perc;
    T = data.AirTemperature1C + 273.15;
    pres = data.AirPressurehPa*100;
    data.RelativeHumidity1Perc = RHice2water(RH_wrt_i,T,pres);

    RH_wrt_i = data.RelativeHumidity2Perc;
    T = data.AirTemperature2C + 273.15;
    pres = data.AirPressurehPa*100;
    data.RelativeHumidity2Perc = RHice2water(RH_wrt_i,T,pres);
    clearvars  RH_wrt_i T pres
end

%% Plotting before removal
    var_name = {'ShortwaveRadiationDownWm2','ShortwaveRadiationUpWm2',...
    'AirTemperature1C','AirTemperature2C','AirTemperature3C','AirTemperature4C',...
    'RelativeHumidity1Perc','RelativeHumidity2Perc','AirPressurehPa',...
    'WindSpeed1ms','WindSpeed2ms','SnowHeight1m','SnowHeight2m'};
    
    var_name_short={'swrd','swru','ta1','ta2','ta3','ta4',...
        'rh1','rh2','ps','ws1','ws2','hs1','hs2'};
    
    f = figure('Visible',vis);%,'outerposition',[0 0 30 40]);
    
    count = 0;
    for i = 1:length(var_name)
        if ismember(var_name{i},data.Properties.VariableNames)
            count = count+1;
        end
    end
        
    ha = tight_subplot(count,1,0.02, [0.05 0.02],[0.07 0.25]);
    
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
            set(f,'CurrentAxes',ha(i)) 

        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('red'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('red'))
        end
        hold on
        if i == 1           
            h_tit = title(station);
            h_tit.Units = 'normalized';
            h_tit.Position(2) = h_tit.Position(2)-2;
            h_tit.Position(1) = h_tit.Position(1)+ 0.6;
        end
        axis tight
        set_monthly_tick(data.time)
        tmp = get(gca,'XTickLabel');
        set(gca,'XTickLabel','')
    end
        set(gca,'XTickLabel',tmp)
    
%% Applying jaws qc code
    ind_qc = find(contains(data.Properties.VariableNames,'qc'));

if ~isempty(ind_qc)
%     qc_name = {data.Properties.VariableNames{ind_qc}}';
%     qc_name(1:4) = [];
    
    qc_name = {'qc_swdn','qc_swup','qc_ttc1','qc_ttc2','qc_tcs1','qc_tcs2',...
        'qc_rh1','qc_rh2','qc_pressure','qc_u1','qc_u2'};

% Assignment of Codes
% �  code = 0 is not used
% �  code = 1 unmodified data
% �  code = 2 linearly interpolated. For snow height missing data replaced
%               with the last available data point.
% �  code = 3 'frozen' e.g. with wind direction when the anemometer is frosted over with the exact
% same value for more than 4 hours.
% �  code = 5 SWin>SWdown
% �  code = 6 corresponds to synthetic wind values. When one WS is not available
%               its value is caluclated from the available level using
%               logarithmic slope-intercept formula is used assuming a 
%               roughness length of 5cm (average condition).
% �  code = 7 corresponds to those cases where aerodynamic theory agreed with the measured logarithmic profile
% above r^2 = 0.97 . In these cases, the theory is selected to predict the 2 and 10 m wind speeds.
% �  code = 8 represents cases when RH data are used to estimate temperatures below -50 C. Note that RH is
% saturated at temperatures less than -45 C and hence the RH values are a direct function of temperature. This
% method has r squared values between 0.8 and 0.98.
% �  code = 9 when temperature values have been corrected for overheating when wind is small and solar
% radiation is great.
% �  code = 9 for some of the 1997 and 1998 Crawford Point snow height data that were synthesized from the
% regression based on the overlap with CP 2 data.

    for i = 1:length(qc_name)       
        data.(var_name{i})(ismember(data.(qc_name{i}),[2])) = NaN;
    end
end
% plotting at this stage
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
            set(f,'CurrentAxes',ha(i)) 

        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('cyan'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('cyan'))
        end
        ylabel(var_name_short{i})
    end
    
%% some filters for unlikely values

data = MaxMinFilter(data, 'AirPressurehPa', 700, 900);

data = MaxMinFilter(data, 'ShortwaveRadiationDownWm2', 0, 950);
data = MaxMinFilter(data, 'ShortwaveRadiationUpWm2', 0, 950);
data.ShortwaveRadiationUpWm2...
    (data.ShortwaveRadiationUpWm2>data.ShortwaveRadiationDownWm2) = NaN;
data.ShortwaveRadiationUpWm2(isnan(data.ShortwaveRadiationDownWm2)) = NaN;

data = MaxMinFilter(data, 'NetRadiationWm2', -500, 10000);
  
T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
data.rh1 = RHice2water(data.RelativeHumidity1Perc,T+273.15,data.AirPressurehPa);
T = nanmean([data.AirTemperature2C, data.AirTemperature4C],2);
data.rh2 = RHice2water(data.RelativeHumidity2Perc,T+273.15,data.AirPressurehPa);
   
data = MaxMinFilter(data, 'rh1', 40, 100);
data = MaxMinFilter(data, 'rh2', 40, 100);
data.RelativeHumidity1Perc(isnan(data.rh1)) = NaN;
data.RelativeHumidity2Perc(isnan(data.rh2)) = NaN;

data = MaxMinFilter(data, 'AirTemperature1C', -70, 40);
data = MaxMinFilter(data, 'AirTemperature2C', -70, 40);
data = MaxMinFilter(data, 'AirTemperature3C', -39, 40);
data = MaxMinFilter(data, 'AirTemperature4C', -39, 40);
data = MaxMinFilter(data, 'AirPressurehPa', 600, 1000);
data = MaxMinFilter(data, 'LongwaveRadiationDownWm2', 50, 400);
data = MaxMinFilter(data, 'WindSpeed1ms', 0, 100);
data = MaxMinFilter(data, 'WindSpeed2ms', 0, 100);

% plotting at this stage
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
            set(f,'CurrentAxes',ha(i)) 

        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('magenta'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('magenta'))
        end
    end
    
%% issue with GCnetwind sensor during a time
    ind = and(data.WindSpeed2ms>9.5,...
    data.time<= datenum('19-Jun-1996'));
    data.WindSpeed2ms(ind)=NaN;
    ind = and(data.WindSpeed1ms>9.5,...
    data.time<= datenum('19-Jun-1996'));
    data.WindSpeed1ms(ind)=NaN;

    % plotting at this stage
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
            set(f,'CurrentAxes',ha(i)) 

        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('light green'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('light green'))
        end
    end
    
%% filter in terms of albedo
data.ShortwaveRadiationUpWm2( ...
    data.ShortwaveRadiationUpWm2 > ...
    0.95 * data.ShortwaveRadiationDownWm2) = NaN;
data.ShortwaveRadiationUpWm2( ...
    data.ShortwaveRadiationUpWm2 < ...
    0.35 * data.ShortwaveRadiationDownWm2) = NaN;

% plotting at this stage
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
            set(f,'CurrentAxes',ha(i)) 

        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('orange'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('orange'))
        end
    end
    
%% periods when WS<0.01ms for more than 6 hours are considered erroneous
ind = data.WindSpeed1ms < 0.01;
no_wind_count = 0;
for i = 1:length(ind)
    if ind(i) == 1
        no_wind_count = no_wind_count +1;
    else
        if no_wind_count>6
            %too long period without wind, leaving flags up
        else
            % gap less than 6 hours putting down the flag
            ind((i-no_wind_count):(i-1)) = 0;
        end
    end
end

data.WindSpeed1ms(ind) = NaN;

ind = data.WindSpeed2ms < 0.01;
no_wind_count = 0;
for i = 1:length(ind)
    if ind(i) == 1
        no_wind_count = no_wind_count +1;
    else
        if no_wind_count>6
            %too long period without wind, leaving flags up
        else
            % gap less than 6 hours putting down the flag
            ind((i-no_wind_count):(i-1)) = 0;
        end
    end
end

data.WindSpeed2ms(ind) = NaN;

% plotting at this stage
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
            set(f,'CurrentAxes',ha(i)) 

        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('turquoise'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('turquoise'))
        end
    end
    
%% measurements less than 0.5m from the ground are discarded
data.WindSpeed1ms(data.WindSensorHeight1m < 0.5)=NaN;
data.WindSpeed2ms(data.WindSensorHeight2m < 0.5)=NaN;
data.AirTemperature1C(data.WindSensorHeight1m < 0.5)=NaN;
data.AirTemperature3C(data.WindSensorHeight1m < 0.5)=NaN;
data.AirTemperature2C(data.WindSensorHeight2m < 0.5)=NaN;
data.AirTemperature4C(data.WindSensorHeight2m < 0.5)=NaN;
data.RelativeHumidity1Perc(data.WindSensorHeight1m < 0.5)=NaN;
data.RelativeHumidity2Perc(data.WindSensorHeight2m < 0.5)=NaN;

% plotting at this stage
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
            set(f,'CurrentAxes',ha(i)) 

        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('pink'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('pink'))
        end
    end
    
%% measurements at unknown heights are discarded
if ~ismember(station,{'SwissCamp','KAN_U','NEEM','NOAA','Miller','KOB'})
    ind_issue = and(isnan(data.WindSensorHeight1m),~isnan(data.SnowHeight1m));
    
    out     = zeros(size(ind_issue'));
    aa      = [0,ind_issue',0];
    ii      = strfind(aa, [0 1]);
    out(ii) = strfind(aa, [1 0]) - ii;

    if max(out) > 3*24*30
%         g= figure('Visible',vis);
%         subplot(2,1,1)
%         hold on
%         plot(data.time, data.WindSensorHeight1m,'LineWidth',1.5)
%         plot(data.time, data.SnowHeight1m,'LineWidth',1.5)
%         legend('Wind sensor height','Snow height')
%         ylabel('Height (m)')
%         set_monthly_tick(data.time)
%         set(gca,'XTickLabel','')
%         axis tight
%         subplot(2,1,2)
%         hold on
%         plot(data.time, data.WindSensorHeight2m,'LineWidth',1.5)
%         plot(data.time, data.SnowHeight2m,'LineWidth',1.5)
%         plot(data.time, ind_issue*5,'LineWidth',1.5)
%         legend('Wind sensor height','Snow height')
%         ylabel('Height (m)')
%         xlabel('Year')
%         set_monthly_tick(data.time)
%         set(gca,'XTickLabelRotation',0)
%         axis tight
%         print(g, sprintf('%s/Height_reconstruction1',OutputFolder), '-dpng')

disp('Warning: some period have snow height but no sensor height.')
        maintenance = ImportMaintenanceFile(station);

        switch station
            case 'NASA-SE'
                date1 = '01-Jan-1998';
                date2 = '01-Jun-2005';
            case 'CP2'
                date1 = '01-Jan-1997';
                date2 = '01-Jun-2001';
        end

        ind_up = find(and(maintenance.date>=datenum(date1),maintenance.date<datenum(date2)));
        disp(datestr(maintenance.date(ind_up)))
        ind = and(data.time>=datenum(date1),data.time<datenum(date2));
        SurfaceHeight1 = -data.SnowHeight1m(ind);
        SurfaceHeight2 = -data.SnowHeight2m(ind);

        for i = 1: length(ind_up)-1
            ind_period = find(and(data.time(ind)>=maintenance.date(ind_up(i)),...
                min(maintenance.date(ind_up(i+1)),data.time(ind)<datenum(date2))));
            SurfaceHeight1(ind_period) = 2 + SurfaceHeight1(ind_period)...
                - SurfaceHeight1(ind_period(find(~isnan(SurfaceHeight1(ind_period)),1,'first')));
            SurfaceHeight2(ind_period) = 3.2 + SurfaceHeight2(ind_period)...
                - SurfaceHeight2(ind_period(find(~isnan(SurfaceHeight2(ind_period)),1,'first')));
        end
    %     SurfaceHeight1(SurfaceHeight1<0) = NaN;
    %     SurfaceHeight2(SurfaceHeight2<0) = NaN;

%         h = figure('Visible',vis);
%         plot(data.time(ind),SurfaceHeight1)
%         hold on
%         plot(data.time(ind),SurfaceHeight2)
%         datetick('x')
%         print(h, sprintf('%s/Height_reconstruction2',OutputFolder), '-dpng')

        data.WindSensorHeight1m(ind) = SurfaceHeight1;
        data.WindSensorHeight2m(ind) = SurfaceHeight2;
    end
end

%% Removing data for which height are unknown
    data.WindSpeed1ms(isnan(data.WindSensorHeight1m ))=NaN;
    data.WindSpeed2ms(isnan(data.WindSensorHeight2m ))=NaN;
    data.AirTemperature1C(isnan(data.WindSensorHeight1m ))=NaN;
    data.AirTemperature3C(isnan(data.WindSensorHeight1m ))=NaN;
    data.AirTemperature2C(isnan(data.WindSensorHeight2m ))=NaN;
    data.AirTemperature4C(isnan(data.WindSensorHeight2m ))=NaN;
    data.RelativeHumidity1Perc(isnan(data.WindSensorHeight1m ))=NaN;
    data.RelativeHumidity2Perc(isnan(data.WindSensorHeight2m ))=NaN;

% plotting at this stage
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
            set(f,'CurrentAxes',ha(i)) 

        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('jaune'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('jaune'))
        end
    end
    
%% Manual removal of suspicious data
% site-specific correction of many erroneous periods
% that is done by checking the data once plotted further down, see if a
% sensor gives strange result then coming back here to correct or
% remove the data
switch station
    case 'NASA-E'
            DV = datevec(data.time);
            ind= ismember(DV(:,2),2:5);
            data.ShortwaveRadiationDownWm2(ind) = NaN;
end
       data = SetErrorDatatoNAN(data,station,vis);
       
	%plotting after removal
    for i = 1:length(var_name)
        if ~ismember(var_name{i},data.Properties.VariableNames)
            continue
        end
        set(f,'CurrentAxes',ha(i)) 
        if ~isempty(strfind(var_name{i},'Relat'))
            T = nanmean([data.AirTemperature1C, data.AirTemperature3C],2);
            temp = RHice2water(data.(var_name{i}),T+273.15,data.AirPressurehPa);
            plot(data.time,temp,'LineWidth',1.5,'Color',RGB('dark green'))
        else
            plot(data.time,data.(var_name{i}),'LineWidth',1.5,'Color',RGB('dark green'))
        end
        datetick('x','yyyy-mm')
        axis tight
        set_monthly_tick(data.time)
        xlim(data.time([1 end]))
        if i <length(var_name)
            set(gca,'XTickLabel','')
        end
    end
        set(gca,'XTickLabelRotation',0)

        set(f,'CurrentAxes',ha(3)) 
        h_leg = legend('jaws qc','max/min','ws issue',...
            'albedo', 'wind freeze', 'too close to ground',...
            'no height','manual removal', 'final' ,'Location','NorthWest');
        h_leg.Position(1) = h_leg.Position(1) +  0.69;
    print(f,sprintf('%s/DataFiltering_%s',OutputFolder,station),'-dpng')

%% smoothing
% see hampel filter and GCnet documentation for further details
%     data.AirTemperature1C = hampel(data.AirTemperature1C,10,1);
%     data.AirTemperature2C = hampel(data.AirTemperature2C,10,1);
%     data.RelativeHumidity1Perc = hampel(data.RelativeHumidity1Perc,10,1);
%     data.RelativeHumidity2Perc = hampel(data.RelativeHumidity2Perc,10,1);
%     data.RelativeHumidity2mPerc = hampel(data.RelativeHumidity2mPerc,10,1);
%     data.ShortwaveRadiationUpWm2 = hampel(data.ShortwaveRadiationUpWm2,5,1);
%     data.ShortwaveRadiationDownWm2 = hampel(data.ShortwaveRadiationDownWm2,5,1);
%     data.NetRadiationWm2 = hampel(data.NetRadiationWm2,5,1);
%     data.AirPressurehPa = hampel(data.AirPressurehPa,48,1); 

%%  Reconstructing temperatures from the two sensors 
% GCnet station have two temperature sensors for each of the two heights
% level. We decide to use first the data from the Thermo Couple sensor and 
% only when this one fails to use the CS500 sensor
        f = figure('Visible',vis);
        ha = tight_subplot(2,2,[0.15 0], [.15 .03], [0 0.02]);
        set(f,'CurrentAxes',ha(1))
        hold on
        scatter( data.AirTemperature1C, data.AirTemperature3C,'.')
        axis tight square
        box on
        ylimit = get(gca,'YLim');
        plot(ylimit,ylimit,'k');
        xlabel('TC level 1')
        ylabel('CS500 level 1')
        
        set(f,'CurrentAxes',ha(3))
        scatter( data.AirTemperature2C, data.AirTemperature4C,'.')
        hold on
        box on
        ylimit = get(gca,'YLim');
        plot(ylimit,ylimit,'k');
        axis tight square
        xlabel('TC level 2')
        ylabel('CS500 level 2')

% at level 1
    set(f,'CurrentAxes',ha(2));
    plot(data.time,data.AirTemperature3C)    
    hold on
    plot(data.time, data.AirTemperature1C)
    axis tight
    xlim([data.time(1) data.time(end)]);
%     set_monthly_tick(data.time); 
    set(gca,'XTickLabel',[],'XMinorTick','on');
    ylabel('Air temperature at\newline level 1 (degC)','Interpreter','tex')
    legend('from CS500','from TC')

% Second Level      
    set(f,'CurrentAxes',ha(4));
    plot(data.time,data.AirTemperature4C)
    hold on
    plot(data.time, data.AirTemperature2C)
    axis tight
    xlim([data.time(1) data.time(end)]);
%     set_monthly_tick(data.time); 
    ylabel('Air temperature at\newline level 2 (degC)','Interpreter','tex')
    datetick('x','yyyy','keeplimits','keepticks')            
    xlabel('Date')
    title(station)
    legend('from CS500','from TC')
    set(gca,'XMinorTick','on');
    print(f, sprintf('%s/T_reconstruction_%s',OutputFolder,station), '-dpng')
    
    % from now AirTemperature1C, resp.AirTemperature2C, contains a combination of the thermocouple
    % and CS500 instruments level 1, resp. 2.
    data.AirTemperature1C(isnan(data.AirTemperature1C)) = ...
        data.AirTemperature3C(isnan(data.AirTemperature1C));
    data.AirTemperature2C(isnan(data.AirTemperature2C)) = ...
        data.AirTemperature4C(isnan(data.AirTemperature2C));

%% Instrument heights
% gaps in air temperature, RH or WS might be at different period. So the
% gap filling process might end up using, at a specific time step,
% temperature from the main station and relative humidity from the
% secondary. We therefor need to keep track of the individual instrument
% heights.
data.WindSensorHeight1m(data.WindSensorHeight1m<0) = NaN;
data.WindSensorHeight2m(data.WindSensorHeight2m<0) = NaN;
data.TemperatureSensorHeight1m(data.TemperatureSensorHeight1m<0) = NaN;
data.TemperatureSensorHeight2m(data.TemperatureSensorHeight2m<0) = NaN;
data.HumiditySensorHeight1m(data.HumiditySensorHeight1m<0) = NaN;
data.HumiditySensorHeight2m(data.HumiditySensorHeight2m<0) = NaN;

% Wind, temperature or humidity sensors that
ind_nan = isnan(data.WindSensorHeight1m );
data.WindSensorHeight1m = hampel(data.WindSensorHeight1m,24*14,0.01);
data.WindSensorHeight1m (ind_nan)=NaN;

ind_nan = isnan(data.WindSensorHeight2m );
data.WindSensorHeight2m = hampel(data.WindSensorHeight2m,24*14,0.01);
data.WindSensorHeight2m (ind_nan)=NaN;

ind_nan = isnan(data.TemperatureSensorHeight1m );
data.TemperatureSensorHeight1m = hampel(data.TemperatureSensorHeight1m,24*14,0.01);
data.TemperatureSensorHeight1m (ind_nan)=NaN;

ind_nan = isnan(data.TemperatureSensorHeight2m );
data.TemperatureSensorHeight2m = hampel(data.TemperatureSensorHeight2m,24*14,0.01);
data.TemperatureSensorHeight2m (ind_nan)=NaN;

ind_nan = isnan(data.HumiditySensorHeight1m );
data.HumiditySensorHeight1m = hampel(data.HumiditySensorHeight1m,24*14,0.01);
data.HumiditySensorHeight1m (ind_nan)=NaN;

ind_nan = isnan(data.HumiditySensorHeight2m );
data.HumiditySensorHeight2m = hampel(data.HumiditySensorHeight2m,24*14,0.01);
data.HumiditySensorHeight2m (ind_nan)=NaN;

%% plot Sensor Height and compare to the reported heights
if ~ismember(station, {'Miller','NOAA'})
    maintenance = ImportMaintenanceFile(station);
    date_change = maintenance.date;

    f = figure('Visible',vis);
    ha = tight_subplot(3,1,.01, [.2 .01], [0.1 0.05]);
   set(ha(1),'Visible','off');
    set(f,'CurrentAxes',ha(2));
    hold on
    plot(data.time, data.WindSensorHeight1m,'b','LineWidth',1.5)
    plot(data.time, data.TemperatureSensorHeight1m,'LineWidth',1.5,'Color',RGB('dark green'),'LineWidth',1.5)
    scatter([datenum(maintenance.date); datenum(maintenance.date)], ...
        [maintenance.W1beforecm/100; maintenance.W1aftercm/100],...
        80,'b','o','filled', 'MarkerFaceColor', 'b')
    scatter([datenum(maintenance.date); datenum(maintenance.date)], ...
        [maintenance.T1beforecm/100; maintenance.T1aftercm/100],...
        80,[51,153,255]/255,'o','filled', 'MarkerFaceColor',RGB('dark green'))

    axis tight
    ylimit=get(gca,'YLim');
    for i = 1:length(date_change)
        h = line(datenum([date_change(i) date_change(i)]),[ylimit(1), ylimit(2)]);
        h.Color = [96,96,96]/255;
        h.LineWidth = 1;
    end
    legendflex({'Wind sensor height from SR',...
        'Temp sensor height from SR',...
        'Wind sensor height from report',...
        'Temp sensor height from report',...
        'Maintenance'}, 'ref', gcf, ...
                           'anchor', {'n','n'}, ...
                           'buffer',[0 -20], ...
                           'ncol',3, ...
                           'fontsize',15,...
                           'title',station,...
                           'Interpreter','none');
    ylabel_obj = ylabel('Height above the surface (m)','Interpreter','tex');
    ylabel_obj.Units = 'Normalized';
    ylabel_obj.Position(2) = ylabel_obj.Position(2)-0.4;
    xlim([data.time(1),data.time(end)])
    set_monthly_tick(data.time); 
    box on
%     h_title = title(station);
%     h_title.Units = 'normalized';
%     h_title.Position(2) = -0.9;
    handle = title('level 1');
     set(handle,'Units','normalized'); 
     set(handle,'Position',[0.95 0.8],'fontsize',15); 
     
     
     if sum(~isnan(data.WindSensorHeight2m))>10 
        set(gca,'XTickLabels',[]);

        set(f,'CurrentAxes',ha(3));
        hold on
        plot(data.time, data.WindSensorHeight2m,'k','LineWidth',1.5)
        scatter(datenum(maintenance.date),maintenance.W2beforecm/100,80,'b','o','filled', 'MarkerFaceColor', 'b')
        scatter(datenum(maintenance.date),maintenance.W2aftercm/100,80,'b','d','filled', 'MarkerFaceColor', 'b')
        scatter(datenum(maintenance.date),maintenance.T2beforecm/100,80,[51,153,255]/255,'o','filled', 'MarkerFaceColor', [51,153,255]/255)
        scatter(datenum(maintenance.date),maintenance.T2aftercm/100,80,[51,153,255]/255,'d','filled', 'MarkerFaceColor', [51,153,255]/255)

        axis tight
        ylimit=get(gca,'YLim');
        for i =1:length(date_change)
            h = line(datenum([date_change(i) date_change(i)]),[ylimit(1), ylimit(2)]);
            h.Color = [96,96,96]/255;
            h.LineWidth = 1;
        end
        handle = title('level 2');
         set(handle,'Units','normalized'); 
         set(handle,'Position',[0.95 0.8],'fontsize',15); 
         box on
         ylabel('')
        xlim([data.time(1),data.time(end)])
        set_monthly_tick(data.time); 
     else
         set(ha(3),'Visible','off')
     end
         
    orient(f,'landscape')
    print(f, sprintf('%s/height_wind_temp_sensors',OutputFolder),'-dpdf')
end

%% keeping useful fields
varlist = {'time',  ...
    'ShortwaveRadiationDownWm2',...
    'ShortwaveRadiationUpWm2',...
    'LongwaveRadiationDownWm2',...
    'LongwaveRadiationUpWm2',...
    'NetRadiationWm2',...
    'AirTemperature1C',...
    'AirTemperature2C',...
    'AirTemperature3C',...
    'AirTemperature4C',...
    'RelativeHumidity1Perc',...
    'RelativeHumidity2Perc',...
    'WindSpeed1ms',...
    'WindSpeed2ms',...
    'WindDirection1deg',...
    'WindDirection2deg',...
    'AirPressurehPa',...
    'SnowHeight1m',...
    'SnowHeight2m',...
    'IceTemperature1C',...
    'IceTemperature2C',...
    'IceTemperature3C',...
    'IceTemperature4C',...
    'IceTemperature5C',...
    'IceTemperature6C',...
    'IceTemperature7C',...
    'IceTemperature8C',...
    'IceTemperature9C',...
    'IceTemperature10C',...
    'WindSensorHeight1m',...
    'WindSensorHeight2m',...
    'Albedo',...
    'ZenithAngledeg',...
    'SZA_jaws',...
    'az',...
    ...'tilt_direction',...
    ...'tilt_angle',...
    'TemperatureSensorHeight1m',...
    'TemperatureSensorHeight2m',...
    'HumiditySensorHeight1m',...
    'HumiditySensorHeight2m'};

ind = [];
for i = 1:length(varlist)
    ind = [ind find(strcmp(data.Properties.VariableNames,varlist{i}))];
end
data = data(:,ind);

%% before anything we interpolate to fill the small gaps
data = InterpTable(data,6);
% data = ResampleTable(data);
end