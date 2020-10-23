clear all
close all
set(0,'defaultfigurepaperunits','centimeters');
   set(0,'DefaultAxesFontSize',16)
   set(0,'defaultfigurecolor','w');
set(0,'defaultfigureinverthardcopy','off');
set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[20 18]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 25]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 25]-0.5]);
addpath(genpath('.\lib'))
addpath(genpath('Input'),genpath('Output'))

WorkingFolder = 'C:\Data_save\2020 JoG data\Corrected';
mkdir(sprintf('%s/Plots',WorkingFolder));
OutputFolder = sprintf('%s/Plots',WorkingFolder);
list = dir(WorkingFolder);
folder_list = {list.name};
folder_list = {folder_list{3:end}};
folder_list(find(strcmp(folder_list,'Plots')))= [];
folder_list(find(strcmp(folder_list,'Plots_saved')))= [];
folder_list(find(strcmp(folder_list,'distrib')))= [];
folder_list(find(strcmp(folder_list,'AWS data')))= [];
folder_list(find(strcmp(folder_list,'Corrected observed firn temperature')))= [];
folder_list(find(strcmp(folder_list,'data out')))= [];
folder_list(find(strcmp(folder_list,'PlotDatasets.m')))= [];
path_list=cell(1,9);
station=cell(1,9);
 count = 0;
 for i = [1 2 4 5 6 7 3 8 9]
    count = count +1;
    path_list{count} = [WorkingFolder '/' folder_list{i}];
    temp = folder_list{i};
    ind =strfind(temp,'_0');
    station{count} = temp(1:ind-1);
end

station2 = station;
for i =1:length(station)
    station2{i} = strrep(station{i},'-','');
    station2{i} = strrep(station2{i},' ','');
end

vis = 'off';

%% Loading subsurface data
    depth_act = cell(1,length(path_list));
    depth_act_2 = cell(1,length(path_list));
    depth_obs_s = cell(1,length(path_list));
    T_ice_mod = cell(1,length(path_list));
    rhofirn = cell(1,length(path_list));
    compaction = cell(1,length(path_list));
    H_surf = cell(1,length(path_list));

    T_ice_obs = cell(1,length(path_list));
    time_mod = cell(1,length(path_list));
    snowbkt = cell(1,length(path_list));
    Surface_Height = cell(1,length(path_list));
    for i = 1:length(path_list)
        disp(station{i})
        % loading the data from the main run
        [depth_act_2{i},depth_act{i}, T_ice_mod{i}, compaction{i}, ...
            time_mod{i}, snowbkt{i}, rhofirn{i},H_surf{i}, depth_obs_s{i}, ...
            T_ice_obs_org{i}, Surface_Height{i}, ~] ...
            = LoadTempModel(path_list{i},station{i});
    end
   disp('Saving')
   tic
    load(strcat(path_list{i},'/run_param.mat'))
depth_obs_save = depth_obs_s;
    toc
 clearvars snowbkt ind coun temp

%% Output surface netcdf files
WriteData = 0;
if WriteData == 1
    for ii = 1:9
        load(strcat(folder_list{ii},'/run_param.mat'))

        % extract surface variables
        namefile = sprintf('%s/%s_surf-bin-%i.nc',path_list{ii},station{ii},1);

        % extract origins
        c.InputAWSFile = sprintf('./Output/Corrected/AWS data/%s/data_%s_combined_hour.txt',...
            station{ii},station{ii});

        [~, ~, ~, ~, ~,...
        ~, T2, ~, z_T2, ~,~, ...
        ~, RH2, ~, z_RH2, ~, ~, ...
        ~, WS2, ~, z_WS2, ~, ~,...
        SRin, SRout, LRin, LRout, ~, ...
        ~, Surface_Height, Tsurf_obs, data_out, c] = ...
        ExtractAWSData(c);
        data_out.Snowfallmweq_origin = isnan(data_out.SurfaceHeightm)*4;  

    % Changing origin coding from:
    % 0: main station
    % 1: CP2 at 6.5 km
    % 2: Swiss Camp at 97 km
    % 3: KAN_U at 67 km
    % 4: RCM
    % 5: MODIS
    % 9: NOAA tower at 2 km
    % 10: Miller's station at 2 km
    %
    % to:
    % 0: main station
    % 0<X<100: adjusted from secondary AWS located X km away from GC-Net station  
    % 101: adjusted from RACMO2.3p2
    % 102: calculated from SRin and nearest daily MODIS albedo (or MODIS albedo
    % climatology if before 2000)

    [data_out.AirTemperature2C_Origin] = ChangeOrigin(data_out.AirTemperature2C_Origin);
    [data_out.RelativeHumidity2_Origin] = ChangeOrigin(data_out.RelativeHumidity2_Origin);
    [data_out.WindSpeed2ms_Origin] = ChangeOrigin(data_out.WindSpeed2ms_Origin);
    [data_out.ShortwaveRadiationDownWm2_Origin] = ChangeOrigin(data_out.ShortwaveRadiationDownWm2_Origin);
    [data_out.ShortwaveRadiationUpWm2_Origin] = ChangeOrigin(data_out.ShortwaveRadiationUpWm2_Origin);
    [data_out.Snowfallmweq_origin] = ChangeOrigin(data_out.Snowfallmweq_origin);

        data = {ncread(namefile,'SRin'), ...
                data_out.ShortwaveRadiationDownWm2_Origin,...
                ncread(namefile,'SRout'), ...
                data_out.ShortwaveRadiationUpWm2_Origin,...
                ncread(namefile,'LRin'), ...
                ncread(namefile,'LRout_mdl'), ...
               ncread(namefile,'SHF'),...
               ncread(namefile,'LHF'),...
               ncread(namefile,'GF'),...
               ncread(namefile,'meltflux'),...
               H_surf{ii},...
               Surface_Height,...
               ncread(namefile,'snowfall')+ncread(namefile,'sublimation'),...
               ncread(namefile,'meltflux')*c.dt_obs/c.dev/c.L_fus/c.rho_water,...   
               ncread(namefile,'sublimation'),...
               ncread(namefile,'snowfall'),...
               data_out.Snowfallmweq_origin,...
               ncread(namefile,'Tsurf')+c.T_0,...
               ncread(namefile,'snowbkt'),...
               ncread(namefile,'theta_2m'),...
               data_out.AirTemperature2C_Origin,...
               ncread(namefile,'RH_2m_wrtw'),...
               data_out.RelativeHumidity2_Origin,...
               ncread(namefile,'ws_10m'),...
               data_out.WindSpeed2ms_Origin};

        varname =  {'SRin','SRin_origin','SRout','SRout_origin',...
            'LRin' 'LRout'...
            'SHF' 'LHF' 'GF' 'meltflux'...
            'H_surf_mod' 'H_surf_obs' 'SMB' 'melt' 'sublimation'...
            'snowfall' 'snowfall_origin'...
            'Tsurf' 'snowbkt',...
            'Ta_2m','Ta_2m_origin',...
            'RH_2m','RH_2m_origin',...
            'WS_10m','WS_10m_origin'};

        unit =  {'Wm-2' '-' 'Wm-2' '-' 'Wm-2' 'Wm-2' 'Wm-2' 'Wm-2' 'Wm-2' 'Wm-2' ...
            'm' 'm' 'm_weq' 'm_weq' 'm_weq' ...
            'm_weq', '-', 'K'  'm_weq','K', '-','%', '-','ms-1', '-'};

        long_varname =  {'Incoming shortwave radiation' ...
            ['Origin of SRin: 0 = main station;'...
            '0<X<100 = adjusted from AWS located X km away;'...
            '101 = adjusted from RACMO2.3p2;'],...
            'Outgoing shortwave radiation' ...
            ['Origin of SRout: 0 = main station;'...
            '102 = calculated from SRin and nearest daily MODIS albedo'...
            '(or MODIS albedo climatology if before 2000)'],...
            'Downward longwave radiation flux (from RACMO2.3p2)'...
            'Calculated upward longwave radiation flux'...
            'Sensible heat flux' ...
            'Latent heat flux' ...
            'Conductive heat flux (>0 when from subsurface to surface)' ...
            'Excess energy available for melt'...
            'Modelled surface height (origin at the bottom of the model)'...
            'Observed surface height (origin at the bottom of the station mast)'...
            'Surface mass balance' ...
            'Melt'...
            'Net sublimation (sublimation if <0, deposition if>0)'...
            'Snowfall' ...
            ['Origin of snowfall estimate: 0 = main station;'...
            '101 = adjusted from RACMO2.3p2;'],...
            'Surface temperature', ...
            'Surface snow bucket',...
            '2m air temperature',...
            ['Origin of Ta: 0 = main station;'...
            '0<X<100 = adjusted from AWS located X km away;'...
            '101 = adjusted from RACMO2.3p2;'],...
            '2m relative humidity with regards to water',...
            ['Origin of RH: 0 = main station;'...
            '0<X<100 = adjusted from AWS located X km away;'...
            '101 = adjusted from RACMO2.3p2;'],...
            '10 m wind speed',...
            ['Origin of WS: 0 = main station;'...
            '0<X<100 = adjusted from AWS located X km away;'...
            '101 = adjusted from RACMO2.3p2;'],...
            };

        switch ii
            case 2
                namefile_out = sprintf('./Output/Corrected/data out/Dye-2_surface.nc');
            case 9
                namefile_out = sprintf('./Output/Corrected/data out/Tunu-N_surface.nc');
            otherwise
                namefile_out = sprintf('./Output/Corrected/data out/%s_surface.nc',station{ii});
        end
        WriteNC_1D(namefile_out, time_mod{ii}, data,...
            varname, unit, long_varname)

    end
     clear data varname unit long_varname
end

%% Output files for Jacqueline Oehri
WriteData = 0;
if WriteData == 1
    station2 = station;
    station2{2} = 'Dye-2';
    station2{6} = 'South Dome';
    station2{9} = 'Tunu-N';

    for ii = 1:9
        load(strcat(folder_list{ii},'/run_param.mat'))

        % extract surface variables
        namefile = sprintf('%s/%s_surf-bin-%i.nc',path_list{ii},station{ii},1);
        switch ii
            case 2
                namefile_out = sprintf('./Output/Corrected/data out/Dye-2_surface.nc');
            case 9
                namefile_out = sprintf('./Output/Corrected/data out/Tunu-N_surface.nc');
            otherwise
                namefile_out = sprintf('./Output/Corrected/data out/%s_surface.nc',station{ii});
        end

        % extract origins
        c.InputAWSFile = sprintf('./Output/Corrected/AWS data/%s/data_%s_combined_hour.txt',...
            station{ii},station{ii});

        [~, ~, ~, ~, ~,...
        ~, T2, ~, z_T2, ~,~, ...
        ~, RH2, ~, z_RH2, ~, ~, ...
        ~, WS2, ~, z_WS2, ~, ~,...
        SRin, SRout, LRin, LRout, ~, ...
        ~, Surface_Height, Tsurf_obs, data_out, c] = ...
        ExtractAWSData(c);
        data_out.Snowfallmweq_origin = isnan(data_out.SurfaceHeightm)*4;  

        % Changing origin coding from:
        % 0: main station
        % 1: CP2 at 6.5 km
        % 2: Swiss Camp at 97 km
        % 3: KAN_U at 67 km
        % 4: RCM
        % 5: MODIS
        % 9: NOAA tower at 2 km
        % 10: Miller's station at 2 km
        %
        % to:
        % 0: main station
        % 0<X<100: adjusted from secondary AWS located X km away from GC-Net station  
        % 101: adjusted from RACMO2.3p2
        % 102: calculated from SRin and nearest daily MODIS albedo (or MODIS albedo
        % climatology if before 2000)

        [data_out.AirTemperature2C_Origin] = ChangeOrigin(data_out.AirTemperature2C_Origin);
        [data_out.RelativeHumidity2_Origin] = ChangeOrigin(data_out.RelativeHumidity2_Origin);
        [data_out.WindSpeed2ms_Origin] = ChangeOrigin(data_out.WindSpeed2ms_Origin);
        [data_out.ShortwaveRadiationDownWm2_Origin] = ChangeOrigin(data_out.ShortwaveRadiationDownWm2_Origin);
        [data_out.ShortwaveRadiationUpWm2_Origin] = ChangeOrigin(data_out.ShortwaveRadiationUpWm2_Origin);
        [data_out.Snowfallmweq_origin] = ChangeOrigin(data_out.Snowfallmweq_origin);

        M = table();
        M.studyloc_abbr = repmat(pad(station2{ii},15), length(data_out.Snowfallmweq_origin),1);
        M.studyloc_nr = repmat(ii, length(data_out.Snowfallmweq_origin),1);

        M.Timestamp = time_mod{ii} - datenum(1900,1,1);
        M.Rnet = ncread(namefile,'SRin') - ncread(namefile,'SRout') + ...
            ncread(namefile,'LRin') - ncread(namefile,'LRout_mdl');
        M.H = ncread(namefile,'SHF');
        M.LE = ncread(namefile,'LHF');
        M.G = ncread(namefile,'GF');
        M.Tair = T2;
        M.Rh = RH2;
        M.Ws = WS2;
        M.Albedo = ncread(namefile,'SRout')./ncread(namefile,'SRin');
        M.Tsurf = ncread(namefile,'Tsurf')+c.T_0;
        M.Kin = ncread(namefile,'SRin');
        M.Kout = ncread(namefile,'SRout');
        M.Lin = ncread(namefile,'LRin');
        M.Lout = ncread(namefile,'LRout_mdl');
        M.Knet = ncread(namefile,'SRin') - ncread(namefile,'SRout') ;
        M.Lnet = ncread(namefile,'LRin') - ncread(namefile,'LRout_mdl');
        M.precip = ncread(namefile,'snowfall');
        M.Tair_height =  z_T2;
        M.RH_height = z_RH2;
        M.Ws_height = z_WS2;
        M.Tair_2m = ncread(namefile,'theta_2m');
        M.RH_2m = ncread(namefile,'RH_2m_wrtw');
        M.Ws_10m	=ncread(namefile,'ws_10m');
        M.Tair_origin = data_out.AirTemperature2C_Origin;
        M.RH_origin = data_out.RelativeHumidity2_Origin;
        M.Ws_origin = data_out.WindSpeed2ms_Origin;
        M.Kin_origin = data_out.ShortwaveRadiationDownWm2_Origin;
        M.Kout_origin=data_out.ShortwaveRadiationUpWm2_Origin;
        M.precip_origin =data_out.Snowfallmweq_origin;

        namefile_out = sprintf('./Output/Corrected/data out/%s.csv',station2{ii});
        writetable(M, namefile_out,'Delimiter',';');
    end
    clear data varname unit long_varname
end

%% Removal of erroneous data

VarName=cell(1,9);
T_ice_obs_fltrd_aux = T_ice_obs_org;
T_ice_obs_fltrd     = T_ice_obs_org;
for i =1:32
    VarName{i} = sprintf('IceTemperature%i',i);
end
for ii = 1:length(time_mod)
    % manual removal of erroneous period
    data1 = array2table(T_ice_obs_fltrd_aux{ii}',...
        'VariableNames',{VarName{1:size(T_ice_obs_fltrd_aux{ii},1)}});
    data1.time = time_mod{ii};
    data2= SetErrorDatatoNAN(data1, station{ii},vis,...
        'Corrected/Plots','noplot');
    T_ice_obs_fltrd_aux{ii} =...
        table2array(data2(:,1:size(T_ice_obs_fltrd_aux{ii},1)))';
    
    % filetering the rest
    for i =1:size(T_ice_obs_fltrd_aux{ii},1)
%          T_ice_obs_org{ii} (i, movvar(T_ice_obs_org{ii}(i,:),24*7,'omitnan')>0.05) = NaN;
%          T_ice_obs_org{ii} (i, movvar(T_ice_obs_org{ii}(i,:),6,'omitnan')>0.1) = NaN;
        T_ice_obs_fltrd_aux{ii} (i, T_ice_obs_fltrd_aux{ii} (i, :)<-40) = NaN;
        if ii == 7
            T_ice_obs_fltrd_aux{ii} (i, T_ice_obs_fltrd_aux{ii} (i, :)>-1) = NaN;
        end
        T_ice_obs_fltrd_aux{ii} (i, T_ice_obs_fltrd_aux{ii} (i, :)<-40) = NaN;
        T_ice_obs_fltrd_aux{ii} (i, :) = hampel(T_ice_obs_fltrd_aux{ii}(i,:),24*7,0.05);      
        T_ice_obs_fltrd{ii} (i, :) = hampel(T_ice_obs_fltrd_aux{ii}(i,:),6,0.1);      
    end
    
    data3 = array2table(T_ice_obs_fltrd{ii}',...
        'VariableNames',{VarName{1:size(T_ice_obs_fltrd{ii},1)}});
    data3.time = time_mod{ii};
    
    PlotRemoval(data1,data2,...
        VarName, station{ii}, vis,'_manual')
    PlotRemoval(data2,data3,...
        VarName, station{ii}, vis,'_filter')
end
clearvars data1 data2 data3 T_ice_obs_fltrd_aux

%% Correcting thermistor location
station3 = station;
station3{1}= 'CP1';
addpath(genpath('../../Data/AWS/'))
T_ice_obs_cor=T_ice_obs_fltrd;
depth_obs_m = depth_obs_s;
 f=figure('Visible',vis);%('outerposition',[1 -1  25 25]);
 set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,0.02,[0.15 0.02],[0.07 0.1]);
for ii = 1:length(path_list)
    disp(station3{ii})
    [depth_obs_s{ii}, T_ice_obs_cor{ii}] = CalculateThermistorDepth(time_mod{ii},...
        depth_act_2{ii}, H_surf{ii},compaction{ii}, depth_obs_save{ii},...
        T_ice_obs_fltrd{ii}, station3{ii});

    % Plotting
    set(f,'CurrentAxes',ha(ii))
    hold on
    plot(time_mod{ii}, H_surf{ii},'LineWidth',2)
    plot(time_mod{ii}, H_surf{ii}-10,'--')

    for kk = 1:size(depth_obs_s{ii},1)
        depth_obs_m{ii}(kk,:) = depth_obs_s{ii}(kk,:)-H_surf{ii}';
        ind_nonan = ~isnan(depth_obs_m{ii}(kk,:));
        plot(time_mod{ii}(ind_nonan),...
                -depth_obs_m{ii}(kk,ind_nonan),...
                '.','Color',RGB(mod(kk,21)+1),'LineWidth',2);
    end
    
    set_plot_style(ii, station, time_mod, 'Height above initial surface level (m)');
end
    print(f, sprintf('%s/depth_therm',OutputFolder), '-dpdf','-bestfit')
    if strcmp(vis,'off')
        close(f)
    end
    
%% Outputting observed subsurface temperature
% disp('Printing corrected observed firn temperature')
% mkdir('.\Output\Corrected\Corrected observed firn temperature')
% clc
% for ii = 1:length(path_list)
%     filename = ['.\Output\Corrected\Corrected observed firn temperature\' ,...
%         station{ii} '_T_firn_obs.nc'];
%     disp(station{ii})
%     tic
%     WriteNC_2D(filename, time_mod{ii}, depth_obs_s{ii},...
%         T_ice_obs_cor{ii}, 'Depth','m',...
%         'Thermistor depth below the surface', ...
%         'T_firn', 'degC', 'Firn temperature'); 
%     toc
% end

%% Ploting observed temperature
% vis = 'on';
 f=figure('Visible',vis);%('outerposition',[1 -1  25 25]);
 set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,[0.02 0.02],[0.15 0.02],[0.07 0.1]);
     load(strcat(path_list{1},'/run_param.mat'))

     scale = 'model';
%      scale = 'station';
ylimit = 20;
    step = 72;

for ii = 1:length(path_list)
 set(f,'CurrentAxes',ha(ii))

    TT = ones(size(depth_obs_s{ii},1),1) * time_mod{ii}';
    [~, ind_bot] = min(abs(depth_obs_s{ii} - (ylimit+1)));
    ind_lim = min(max(ind_bot)+1, size( T_ice_obs_cor{ii},1));
    switch scale
        case 'model'
            depth_temp = depth_obs_m{ii};

        case 'station'
            depth_temp = depth_obs_s{ii};
            for i = 1: size(depth_temp,1)
                depth_temp(i,:) = depth_temp(i,:) - Surface_Height{ii}';
            end
    end

    col = PlotTemp(TT(1:ind_lim, 1:step:end),...
        depth_temp(1:ind_lim, 1:step:end),...
        T_ice_obs_cor{ii}(1:ind_lim, 1:step:end),...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','jet',...
        'Interp','on',...
        'XLabel','',...
        'YLabel',' ',...
        'CLabel','Firn temperature (^oC)',...
        'Range', -40:1:0,...
        'FlatSurface','no');
%     col.FontSize = 12;
    if ii == 1
        col.Position(1) = 0.92;
        col.Position(4) = 0.9 ;
        col.Position(2) = 0.08;
    else
        col.Position(1) = col.Position(1)+ 3;
    end
    
    switch scale
        case 'model'
            plot(time_mod{ii}, - H_surf{ii},'LineWidth',2)
            plot(time_mod{ii}, 10 - H_surf{ii},':','LineWidth',2)
        case 'station'
            plot(time_mod{ii}, - Surface_Height{ii}','LineWidth',2)
            plot(time_mod{ii}, 10 - Surface_Height{ii}',':','LineWidth',2)
    end
    
%     col.TickLength=col.TickLength*5;
    temp = col.YTickLabel;
    for k = 1:length(temp)
        if (k-1)/5==floor((k-1)/5)
        else
            temp(k,:)=' ';
        end
    end
    col.YTickLabel=temp;
    
    set_plot_style(ii, station, time_mod, 'Height above initial surface level (m)');
    ylim([- H_surf{ii}(end) - 0.5 15 + 0.5])

 ha(ii).Position(3) = 0.25;
 ha(ii).Position(4) = 0.27;
 %    set (gca,'XTick',1995:2018,'XMinorTick','off','YMinorTick','on','XTickLabelRotation',45)
end

print(f, sprintf('%s/T_ice_obs_cor',OutputFolder), '-dpdf','-bestfit')
    if strcmp(vis,'off')
        close(f)
    end
 clearvars TT ha temp col
   
%% Ploting modelled temperature
 f=figure('Visible',vis);%('outerposition',[1 -1  25 25]);
 set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,0.02,[0.15 0.02],[0.07 0.1]);
     load(strcat(path_list{1},'/run_param.mat'))
ylimit = 40;
for ii = 1:length(path_list)
 set(f,'CurrentAxes',ha(ii))

    TT = ones(size(depth_act{ii},1),1) * time_mod{ii}';
    [~, ind_bot] = min(abs(depth_act{ii} - (ylimit+1)));
    ind_lim = min(max(ind_bot)+1, size( T_ice_mod{ii},1));
    step = 72;
    
    col = PlotTemp([TT(1, 1:step:end); TT(1:ind_lim, 1:step:end)],...
        [-H_surf{ii}(1:step:end)'; depth_act{ii}(1:ind_lim, 1:step:end)],...
        [T_ice_mod{ii}(1, 1:step:end); T_ice_mod{ii}(1:ind_lim, 1:step:end)]-273.15,...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','jet',...
        'Interp','on',...
        'XLabel','',...
        'YLabel',' ',...
        'CLabel','Firn temperature (^oC)',...
        'Range', -40:1:0,...
        'FlatSurface','no');
%     col.FontSize = 12;
    if ii == 1
        col.Position(1) = 0.92;
        col.Position(4) = 0.9 ;
        col.Position(2) = 0.08;
    else
        col.Position(1) = col.Position(1)+ 3;
    end
    plot(time_mod{ii},-H_surf{ii},'LineWidth',2)
    plot(time_mod{ii},10-H_surf{ii},':','LineWidth',2)
    
%     col.TickLength=col.TickLength*5;
    temp = col.YTickLabel;
    for k = 1:length(temp)
        if (k-1)/5==floor((k-1)/5)
        else
            temp(k,:)=' ';
        end
    end
    col.YTickLabel=temp;
    
    set_plot_style(ii, station, time_mod, 'Height above initial surface level (m)');
%     ylim([- H_surf{ii}(end) - 0.5 max(max(depth_obs_m{ii}))+0.5])
    ylim([- H_surf{ii}(end) - 0.5 15 + 0.5])

 ha(ii).Position(3) = 0.25;
 ha(ii).Position(4) = 0.27;
 %    set (gca,'XTick',1995:2018,'XMinorTick','off','YMinorTick','on','XTickLabelRotation',45)
end

	    print(f, sprintf('%s/T_ice_mod',OutputFolder), '-dpdf','-bestfit','-bestfit')
% save(sprintf('%s/result_normal',OutputFolder),'depth_act','T_ice_mod','-v7.3')

    if strcmp(vis,'off')
        close(f)
    end
clearvars TT ha temp col
    
%% Ploting modelled density
 f=figure('Visible',vis);%('outerposition',[1 -1  25 25]);
 set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,0.02,[0.15 0.02],[0.07 0.1]);
     load(strcat(path_list{1},'/run_param.mat'))
ylimit = 40;
for ii = 1:length(path_list)
 set(f,'CurrentAxes',ha(ii))

    TT = ones(size(depth_act{ii},1),1) * time_mod{ii}';
    [~, ind_bot] = min(abs(depth_act{ii} - (ylimit+1)));
    ind_lim = min(max(ind_bot)+1, size( T_ice_mod{ii},1));
    step = 72;
    
    col = PlotTemp([TT(1, 1:step:end); TT(1:ind_lim, 1:step:end)],...
        [-H_surf{ii}(1:step:end)'; depth_act{ii}(1:ind_lim, 1:step:end)],...
        [rhofirn{ii}(1, 1:step:end); rhofirn{ii}(1:ind_lim, 1:step:end)],...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','hsv',...
        'Interp','on',...
        'XLabel','',...
        'YLabel',' ',...
        'CLabel','Firn density (kg m-3)',...
        'Range', 315:10:900,...
        'FlatSurface','no');
%     col.FontSize = 12;
    if ii == 1
        col.Position(1) = 0.92;
        col.Position(4) = 0.9 ;
        col.Position(2) = 0.08;
    else
        col.Position(1) = col.Position(1)+ 3;
    end
    plot(time_mod{ii},-H_surf{ii},'LineWidth',2)
    plot(time_mod{ii},10-H_surf{ii},':','LineWidth',2)
    
%     col.TickLength=col.TickLength*5;
    temp = col.YTickLabel;
    for k = 1:length(temp)
        if (k-1)/5==floor((k-1)/5)
        else
            temp(k,:)=' ';
        end
    end
    col.YTickLabel=temp;
    
    set_plot_style(ii, station, time_mod, 'Height above initial surface level (m)');
%     ylim([- H_surf{ii}(end) - 0.5 max(max(depth_obs_m{ii}))+0.5])
    ylim([- H_surf{ii}(end) - 0.5 20+0.5])

 ha(ii).Position(3) = 0.25;
 ha(ii).Position(4) = 0.27;
 %    set (gca,'XTick',1995:2018,'XMinorTick','off','YMinorTick','on','XTickLabelRotation',45)
end

	    print(f, sprintf('%s/rho_mod',OutputFolder), '-dpdf','-bestfit')
% save(sprintf('%s/result',OutputFolder),'depth_act','T_ice_mod','-v7.3')
    if strcmp(vis,'off')
        close(f)
    end
    clearvars TT ha temp col

%% Interpolating modelled firn temperature at observation depth
T_subsurf_mod = cell(1,9);
for ii = 1:length(path_list)
       disp(station{ii})
tic
    T_subsurf_mod{ii} = NaN(size(depth_obs_s{ii}));
    %interpolating model at observed depth
    for j = 1:length(time_mod{ii})
        T_subsurf_mod{ii}(:,j) = interp1( depth_act_2{ii}(:,j),...
            T_ice_mod{ii}(:,j),  ...
            depth_obs_s{ii}(:,j))-c.T_0;
    end
    ind_in_first_layer = depth_obs_s{ii} ...
        <= depth_act_2{ii}(1,:);
    T_subsurf_mod{ii}(ind_in_first_layer) = NaN;
        plot(time_mod{ii},T_subsurf_mod{ii}(1,:),'LineWidth',2)
  toc
end

%% Ploting temperature difference

 f=figure('Visible',vis);%('outerposition',[1 -1  25 25]);
 set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,0.02,[0.15 0.02],[0.07 0.1]);
     load(strcat(path_list{1},'/run_param.mat'))

ylimit = 20;
for ii = 1:length(path_list)
 set(f,'CurrentAxes',ha(ii))

    TT = ones(size(depth_obs_s{ii},1),1) * time_mod{ii}';
    [~, ind_bot] = min(abs(depth_obs_s{ii} - (ylimit+1)));
    ind_lim = min(max(ind_bot)+1, size( T_ice_obs_cor{ii},1));
    step = 24;
    
% Subsurface temperature bias
T_diff = T_subsurf_mod{ii} - T_ice_obs_cor{ii};

    col = PlotTemp(TT(1:ind_lim, 1:step:end),...
        depth_obs_m{ii}(1:ind_lim, 1:step:end),...
        T_diff(1:ind_lim, 1:step:end),...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','bwr_cmap',...
        'Interp','on',...
        'XLabel','',...
        'YLabel',' ',...
        'CLabel','Modelled minus observed firn temperature (^oC)',...
        'Range', -5:0.5:5,...
        'FlatSurface','no');
%     col.FontSize = 12;
    if ii == 1
        col.Position(1) = 0.92;
        col.Position(4) = 0.9 ;
        col.Position(2) = 0.08;
    else
        col.Position(1) = col.Position(1)+ 3;
    end
    plot(time_mod{ii}, - H_surf{ii},'LineWidth',2)
    
    temp = col.YTickLabel;
    for k = 1:length(temp)
        if (k-1)/5==floor((k-1)/5)
        else
            temp(k,:)=' ';
        end
    end
    col.YTickLabel=temp;

    set_plot_style(ii, station, time_mod, 'Height above initial surface level (m)');
    ylim([- H_surf{ii}(end) - 0.5 15 + 0.5])

 ha(ii).Position(3) = 0.25;
 ha(ii).Position(4) = 0.27;
end
	    print(f, sprintf('%s/T_ice_diff',OutputFolder), '-dpdf','-bestfit')
    if strcmp(vis,'off')
        close(f)
    end
clearvars TT ha temp col

%% Calculating observed and modelled 10 m firn temperature
depth_deep = 10;
disp('Calculating observed and modelled 10 m firn temperature')
tic
T_deep_obs = cell(1,9);

for ii = 1:length(path_list)
    disp(ii)
    T_deep_obs{ii} = zeros(size(time_mod{ii}));

    for i =1:length(time_mod{ii})        
        ind = and( and(~isnan(depth_obs_s{ii}(:,i)), ...
             ~isnan(T_ice_obs_cor{ii}(:,i))), ...
            depth_obs_s{ii}(:,i)~=0);
        if sum(ind)>3 && max(depth_obs_s{ii}(ind,i))>5
            [~, index, ~] = unique(depth_obs_s{ii}(ind,i),'last');
            indexToDupes = ...
                find(not(ismember(1:numel(depth_obs_s{ii}(ind,i)),index)));
            if ~isempty(indexToDupes)
                ind(find(ind) == indexToDupes) = 0;
            end
            if max(depth_obs_s{ii}(ind(~isnan(T_ice_obs_cor{ii}(ind,i))),i)) ...
                    > depth_deep - 1.5
                T_deep_obs{ii}(i) = interp1(depth_obs_s{ii}(ind,i),...
                    T_ice_obs_cor{ii}(ind,i),depth_deep,'linear','extrap');
            else
                T_deep_obs{ii}(i) = interp1(depth_obs_s{ii}(ind,i),...
                    T_ice_obs_cor{ii}(ind,i),depth_deep,'linear');                
            end
        else
            T_deep_obs{ii}(i)=NaN;
        end
    end
end

%% modelled temperature at 10 m
T_deep_mod = cell(1,9);

for ii = 1:length(path_list)
    disp(ii)
    T_deep_mod{ii} = zeros(1,size(depth_act_2{ii},2));
     [depth_alt, M_alt, avg_shallow] = ...
        InterpGridColumnWise(depth_act_2{ii}, T_ice_mod{ii}, 10,...
        'lin', 'VolumeWeighted');
    ind_mat = depth_alt==10;
    ind_double = find(sum(ind_mat,1)>1);
    for k = 1:length(ind_double);
        ind_first_double = find(ind_mat(: , ind_double(k)) ,1,'first');
        ind_mat( ind_first_double , ind_double(k)) = 0;
    end
    T_deep_mod{ii} = M_alt(ind_mat);
end
toc

%% Fig. 5: Plotting observed vs. modelled 10 m firn temperature

f=figure('Visible',vis,'Renderer','Painter','PaperSize',[22 16]);%('outerposition',[1 -1  25 25]);
set(gcf,'Position',[0    1.0583   22   13])
 ha = tight_subplot(3,3,[0.03 0.01],[0.15 0.04],[0.12 0.05]);

for ii = 1:length(path_list)
    % plotting
    set(f,'CurrentAxes',ha(ii))
    hold on
    plot(time_mod{ii},interp1gap(T_deep_obs{ii},24),'LineWidth',3,'Color',RGB('dark gray'))
%     plot(time_mod{ii},T_deep_obs{ii},'LineWidth',3,'Color',RGB('dark gray'))
    
    plot(time_mod{ii},T_deep_mod{ii}-c.T_0,'LineWidth',2,'Color',RGB('black'))
    lm = fitlm(time_mod{ii},T_deep_mod{ii}-c.T_0);
    plot(time_mod{ii}, ...
        lm.Coefficients.Estimate(1)+lm.Coefficients.Estimate(2)*time_mod{ii},...
        'Color',RGB('red'),...
        'LineWidth',2)
    disp(station{ii})
    disp(lm.Coefficients.Estimate(2)*365*10)

    set_plot_style(ii, station, time_mod,...
        '10 m firn temperature (^{o}C)');
    switch ii
        case {1,2,3} 
            ylim([-24 -12])
            a = -22;
            set(gca,'YTick',-25:5:-10)
        case {4,5,6}
            ylim([-26 -17])
            a = -25;
            set(gca,'YTick',-25:5:-10)
        case {7,8,9}
            ylim([-34 -26])
            a = -33;
    end
    if ismember(ii,[7 9])
        xlabel('')
    end
    text(datenum(2000.5,1,1),a,...
        sprintf('trend: %0.2f ^oC decade^{-1}',lm.Coefficients.Estimate(2)*365*10),...
        'FontWeight','bold','Interpreter','tex',...
        'Color',RGB('red'),...
        'FontSize',12)
    if ismember(ii,[2 3 5 6 8 9])
        set(gca,'YTickLabel','')
    end
    ha(ii).TickLength = ha(ii).TickLength*2;
    xlim(datenum([1998 2018],1,1))
%     set(gca,'XMinorTick','off')
end
legendflex({'Observed','Modelled','Linear trend in modelled firn temperature'},...
    'nrow',1,...
    'ref',gcf,...
    'anchor',{'n','n'},...
    'box','off',...
    'buffer',[0 8]);

set(gcf,'PaperPositionMode','auto');         
% set(gcf,'PaperOrientation','landscape');
print(f, sprintf('%s/JOG-19-0112.Figure5.pdf',OutputFolder), '-dpdf','-bestfit')
    if strcmp(vis,'off')
        close(f)
    end
    
%% Scatter plot at each station

summary_fit_therm = array2table(NaN(3,length(path_list)));
summary_fit_therm.Properties.VariableNames = station2;
summary_fit_therm.Properties.RowNames = {'R2_T_ice' 'RMSE_T_ice' 'ME_T_ice'} ;

 f=figure('Visible',vis,'outerposition',[1 -1  20 20]);
%  set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,[0.03 0.01],[0.15 0.02],[0.08 0.11]);
for ii = 1:length(path_list)
set(f,'CurrentAxes',ha(ii))
disp(station{ii})
hold on
    for kk = 1:size(T_ice_obs_cor{ii},1)
         scatter(T_ice_obs_cor{ii}(kk,:), T_subsurf_mod{ii}(kk,:), 5,'fill')
    end

    axis tight
    xlimit = get(gca, 'Xlim');
    ylimit = get(gca, 'Ylim');
    plot([-40 0],[-40 0],'--k','LineWidth',2)
    
    lm = fitlm(T_ice_obs_cor{ii}(:),T_subsurf_mod{ii}(:));
    clearvars R2_T_ice RMSE_T_ice ME_T_ice
    summary_fit_therm.(station2{ii})(1) = lm.Rsquared.Adjusted;
    summary_fit_therm.(station2{ii})(2) = lm.RMSE;
    summary_fit_therm.(station2{ii})(3) = nanmean(T_subsurf_mod{ii}(:)-T_ice_obs_cor{ii}(:));

    ylim([-40 0])
    xlim([-40 0])
    xlimit = get(gca, 'Xlim');
    ylimit = get(gca, 'Ylim');
    y = ylimit(2)- abs(ylimit(2)-ylimit(1))/7;
    x = xlimit(1)+ abs(xlimit(2)-xlimit(1))/10;
    text(x,y,sprintf('R^2 = %0.02f\nRMSE = %0.02f\nME = %0.02f',...
         summary_fit_therm.(station2{ii})(1), ...
         summary_fit_therm.(station2{ii})(2) ,...
         summary_fit_therm.(station2{ii})(3)),'Interpreter','tex');
     
     axis square
    box on
    set(gca,'XMinorTick','on','YMinorTick','on');

    h_tit = title([char(96+ii) ' ' station3{ii}]);
    h_tit.FontSize = 12;
    h_tit.Units = 'Normalized';
    h_tit.Position(2) = 0.05;
    if ii == 8
        xlabel('Observed firn temperature (^oC)','Interpreter','tex')
    end
    if ii==4
        ylabel('Modelled firn temperature (^oC)','Interpreter','tex')
    end
    if ismember(ii,1:6)
        set(gca,'XTickLabel','')
    end
    if ismember(ii,[2 3 5 6 8 9])
        set(gca,'YTickLabel','')
    end
    ha(ii).TickLength = 2*ha(ii).TickLength;
grid on
end
legendflex({'#1',...
    '#2',...
    '#3',...
    '#4',...
    '#5',...
    '#6',...
    '#7',...
    '#8',...
    '#9',...
    '#10'},...
    'ref',gcf,...
    'anchor',{'e','e'},...
    'box','off',...
    'title','Sensor:')
    
writetable2csv(summary_fit_therm,sprintf('%s/summary_fit_therm.txt',OutputFolder))

 print(f, sprintf('%s/T_firn_validation',OutputFolder), '-dpdf','-bestfit')
if strcmp(vis,'off')
    close(f)
end

%% Ploting comparison for each thermistor
for ii = 1:length(path_list)

    f = figure('Visible',vis);
     ha = tight_subplot(10,1,0.02,[0.08 0.05],[0.07 0.1]);

    for kk = 1:10
        set(f,'CurrentAxes',ha(kk))
        hold on
        T_subsurf_mod{ii}(kk,isnan(T_ice_obs_cor{ii}(kk,:)))=NaN;
        plot(time_mod{ii}, T_ice_obs_cor{ii}(kk,:))
        plot(time_mod{ii}, T_subsurf_mod{ii}(kk,:))
        h_text = text(time_mod{ii}(24*60), ...
           -18,...
            sprintf('%i',kk));
        h_text.FontSize = 15;
        h_text.FontWeight = 'bold';
        h_text.Units = 'Normalized';
        h_text.Position(1:2) = [0.02 0.5];
        set_monthly_tick(time_mod{ii})
        set(gca,'XTickLabelRotation',0)
        axis tight
        xlim(time_mod{ii}([1 end]))
        if ismember(kk,1:9)
           set(gca,'XTickLabel','')
        else
           xlabel('Year')
        end
        if kk == 1
            h_tit = title(station{ii});
            h_tit.Units = 'normalized';
            h_tit.Position(2) = h_tit.Position(2)-0.5;
            legendflex({'Observed','Modelled'},...
                'ref',gcf,'anchor',{'n' 'n'},'nrow',1)
        end
        if ismember(kk,2:2:size(T_ice_obs_cor{ii},1))
           set(gca,'YAxisLocation','right')
        end
        if kk==5
            ylabel('Firn temperature (degC)')
        end
    end

    print(f, sprintf('%s/T_firn_validation2_%s',OutputFolder,station{ii}), '-dpdf','-bestfit')
    if strcmp(vis,'off')
        close(f)
    end

end    

%% Top 20 m cold content
ColdContent = @(rho,thickness,T_firn) ...
    2.108 .* rho .* thickness .* (c.T_0 - T_firn);
%    kJ/kg/K * kg/m3 * m            * K
% in kJ/m2

 CC = cell(1,9);
 CC_20 = cell(1,9);
 T_avg_20 = cell(1,9);
 rho_avg_20 = cell(1,9);
summary_CC = NaN(9,9);

 for ii = 1:length(path_list)
    disp(station{ii})
    tic
    [depth_alt, rho_alt, rho_20] = ...
        InterpGridColumnWise(depth_act_2{ii}, rhofirn{ii}, 20, 'next','VolumeWeighted',[]);
    [~, T_ice_alt, T_20] = ...
        InterpGridColumnWise(depth_act_2{ii}, T_ice_mod{ii}, 20, 'lin','MassWeighted', rho_alt);
   

    th  = depth_act_2{ii};
    th(2:end,:)  = depth_act_2{ii}(2:end,:)-depth_act_2{ii}(1:end-1,:);
        
    thickness_alt  = depth_alt;
    thickness_alt(2:end,:)  = depth_alt(2:end,:)-depth_alt(1:end-1,:);
    

    CC{ii} = ColdContent(rhofirn{ii},th,T_ice_mod{ii});     
    CC_20{ii}  = ColdContent(rho_20, 20*ones(size(rho_20)), T_20);
    T_avg_20{ii} = T_20;
    rho_avg_20{ii} = rho_20;
    
    lm = fitlm(time_mod{ii}, CC_20{ii});
    summary_CC(1,ii) = mean(CC_20{ii});
    summary_CC(2,ii) = lm.Coefficients.Estimate(2)*10;
    summary_CC(3,ii) = coefTest(lm);
    lm = fitlm(time_mod{ii}, T_avg_20{ii}-273.15);
    summary_CC(4,ii) = mean(T_avg_20{ii}-273.15);
    summary_CC(5,ii) = lm.Coefficients.Estimate(2)*10;
    summary_CC(6,ii) = coefTest(lm);
    lm = fitlm(time_mod{ii}, rho_avg_20{ii});
    summary_CC(7,ii) = mean(rho_avg_20{ii});
    summary_CC(8,ii) = lm.Coefficients.Estimate(2)*10;
    summary_CC(9,ii) = coefTest(lm);

    toc
end
close all
dlmwrite( sprintf('%s/summary_CC_20.csv',OutputFolder), summary_CC,';')

%% Plotting CC20
f=figure('Visible',vis,'outerposition',[1 -1  23 15]);
%  ha = tight_subplot(3,3,[0.02 0.05],[0.15 0.02],[0.07 0.1]);
 ha = tight_subplot(1,2,0.02,0.2,[0.1 0.18]);
ps = @(m,rho) max(0,m.*(1./rho -1/917));

FAC20 = cell(1,9);
FC20 = cell(1,9);
col = linspecer(9);
for ii = 1:length(path_list)
   
    set(f,'CurrentAxes',ha(1))
    hold on
    plot(time_mod{ii}, CC_20{ii}/334,'LineWidth',2,'Color',col(ii,:))

        set(f,'CurrentAxes',ha(1))

    FAC20{ii} = ps(20*rho_avg_20{ii}, rho_avg_20{ii});
    FC20{ii} = max(0, 20*843 - 917*(20-FAC20{ii}));

        set(f,'CurrentAxes',ha(2))
    hold on
    plot(time_mod{ii}, FC20{ii}/ FC20{ii}(1)*100,'LineWidth',2,'Color',col(ii,:))
    disp(station{ii})
    disp((FC20{ii}(end)-FC20{ii}(1))/FC20{ii}(1))
end
    set(f,'CurrentAxes',ha(1))
axis tight
ylabel('Refreezing capacity in top 20 m (% relative to June 1st, 1998))')
xlabel('Year')
box on
ha(1).XTick = datenum(1994:2:2018,1,1);
ha(1).XMinorTick = 'on';
datetick('x','yyyy','keepticks','keeplimits')
ha(1).XTickLabelRotation = 45;

    set(f,'CurrentAxes',ha(2))
axis tight
ylabel('Firn capacity in top 20 m (mm)','Interpreter','tex')
box on
xlabel('Year')
ha(2).XTick = datenum(1994:2:2018,1,1);
ha(2).XMinorTick = 'on';
ha(2).XTickLabelRotation = 45;
datetick('x','yyyy','keepticks','keeplimits')
ha(2).YAxisLocation = 'right';
    ha(1).TickLength = 2*ha(1).TickLength;
    ha(2).TickLength = 2*ha(2).TickLength;

legendflex({station{[1:6]}},'ref',gcf,'anchor',{'n','n'},'nrow',3)

print(f, sprintf('%s/ColdContent_20_all',OutputFolder), '-dpdf','-bestfit')


% %% Plotting T20
% f=figure('Visible',vis,'outerposition',[1 -1  25 25]);
%  ha = tight_subplot(3,3,[0.02 0.05],[0.15 0.02],[0.07 0.1]);
% %  ha = tight_subplot(1,1,0.02,[0.15 0.02],[0.07 0.1]);
% 
% for ii = 1:length(path_list)
%    
%     set(f,'CurrentAxes',ha(ii))
% hold on
%     plot(time_mod{ii}, T_avg_20{ii}-273.15,'LineWidth',2)
%     set_plot_style(ii, station, time_mod, '20 m average firn temperature (^oC)');
% axis fill
%         ylim(floor([min(T_avg_20{ii})-2 max(T_avg_20{ii})+2])-273.15)
% end
% print(f, sprintf('%s/T_20',OutputFolder), '-dpdf','-bestfit')

%% Density validation

if exist('../matlab_functions/Core_all.mat', 'file') == 2
    disp('Plotting comparison with density profiles from cores')
    
    load ../matlab_functions/Core_all.mat
    
    avg_dens_mod = [];
    avg_dens_meas = [];
for ii=1:length(station)
    i_core = FindCore(Core,'NearestCodeLocation',station{ii});
%     CoreList({Core{i_core}})
    ind_remove = zeros(size(i_core));
    for j = 1:length(i_core)
        if Core{i_core(j)}.Info.DateCored.Year <1999 ||...
                Core{i_core(j)}.Info.DateCored.Year > 2018
            ind_remove(j) = 1;
        end
        if strcmp(Core{i_core(j)}.Info.Name,'core_9_2016')
            ind_remove(j) = 1;
        end 
    end
    i_core(ind_remove==1) = [];
    
    % ordering cores in chronological order
    dates = zeros(size(i_core));
    disp(station{ii})
    for i = 1:length(i_core)
        disp([Core{i_core(i)}.Info.Name '    '])
        dates(i) = datenum(Core{i_core(i)}.Info.DateCored);
    end
    
    [~, i_ordered] = sort(dates);
    i_core = i_core(i_ordered);
    
    thickness_act = depth_act_2{ii};
    thickness_act(2:end,:) = depth_act_2{ii}(2:end,:) - depth_act_2{ii}(1:end-1,:);

    num_plot=8;
    ylim_core = 30;
    
    density_meas{ii} = NaN(4,max(1,length(i_core)));
    density_mod{ii} = NaN(4,max(1,length(i_core)));

    if ~isempty(i_core)
        f = figure('Visible','off');
        [ha, ~] = tight_subplot(1, num_plot, 0.01, [0.12 0.01], [0.05 0.01]);
        count = 0;

        for jj = i_core
            count = count+1;
            if count <= num_plot
                set(f,'CurrentAxes',ha(count))
            else
                i_file = 1;
                NameFile = sprintf('%s/CorePlot_%s_%i.tif',OutputFolder,station{ii}, i_file);
                while exist(NameFile, 'file') == 2
                    i_file = i_file + 1;
                NameFile = sprintf('%s/CorePlot_%s_%i.tif',OutputFolder,station{ii}, i_file);
                end
                print(f,NameFile,'-dpdf');
                if strcmp(vis,'off')
                    close(f);
                end
                f = figure('Visible','off');
                [ha, ~] = tight_subplot(1, num_plot, 0.01, [0.12 0.01], [0.05 0.01]);
                count = 1;
                set(f,'CurrentAxes',ha(count))
            end

            time_core = datenum(Core{jj}.Info.DateCored);

            temp = abs(time_mod{ii} - time_core);
            [~, ind_time] = min(temp);
            Core_2{1} = Core{jj};

            for jk = 2:length(depth_act_2{ii}(:,ind_time))
                if thickness_act(jk,ind_time)>0.20
                    ind_depth = find(...
                        and(Core_2{1}.Data.Depth < depth_act_2{ii}(jk,ind_time)*100,...
                        Core_2{1}.Data.Depth > depth_act_2{ii}(jk-1,ind_time)*100));
                    Core_2{1}.Data.Density(ind_depth) = nanmean(Core_2{1}.Data.Density(ind_depth));
                    if ~isempty(Core_2{1}.Data.Type_perc)
                        Core_2{1}.Data.Type_perc...
                            (min(length(Core_2{1}.Data.Type_perc),ind_depth)) ...
                            = nanmean(Core_2{1}.Data.Type_perc(min(length(Core_2{1}.Data.Type_perc),ind_depth)));
                    end
                end
            end

            Core_2{2}.Info.Densities = 'y';
            Core_2{2}.Info.DateCored = datetime(datestr(time_mod{ii}(ind_time)));
            Core_2{2}.Info.Name = 'Model';
            Core_2{2}.Info.NearestCodeLocation = station{ii};

            depth_max = floor(depth_act_2{ii}(end,ind_time)*100);
            Core_2{2}.Data.Density = zeros(depth_max,1);
            Core_2{2}.Data.Depth  = [1:depth_max]';
            Core_2{2}.Data.Type = cell(depth_max,1);
            Core_2{2}.Data.Type_perc = zeros(depth_max,1);

            for i = length(depth_act_2{ii}(:,ind_time)):-1:1
                ind = find(Core_2{2}.Data.Depth<=100*depth_act_2{ii}(i,ind_time));
                Core_2{2}.Data.Density(ind) = rhofirn{ii}(i,ind_time);
                for j = 1:length(ind)
                    Core_2{2}.Data.Type{ind(j)} = 'firn';
                end
    %             Core_2{2}.Data.Type_perc(ind) = snic(i,ind_time)/thickness_weq(i,ind_time)*100;
            end

            ind_up = find(and(~isnan(Core{jj}.Data.Density), ...
                Core{jj}.Data.Depth<=500));
            
            ind_mid = find(and(~isnan(Core{jj}.Data.Density),...
                and(Core{jj}.Data.Depth>500,Core{jj}.Data.Depth<2000)));
            % case where the model is shorter than the core
            ind = ind_mid>length(Core_2{2}.Data.Density);
            ind_mid(ind)=[]; 
            
            ind_down = find(and(~isnan(Core{jj}.Data.Density),...
                Core{jj}.Data.Depth>=2000));
            % case where the model is shorter than the core
            ind = ind_down>length(Core_2{2}.Data.Density);
            ind_down(ind)=[];
            
            ind_nonan=find(~isnan(Core{jj}.Data.Density));
            ind_nonan(ind_nonan>length(Core_2{2}.Data.Density))=[];

            density_meas{ii}(:,jj) = [mean(Core{jj}.Data.Density(ind_up));...
                mean(Core{jj}.Data.Density(ind_mid));...
                mean(Core{jj}.Data.Density(ind_down));...
                mean(Core{jj}.Data.Density(ind_nonan))];

            density_mod{ii}(:,jj) = [mean(Core_2{2}.Data.Density(ind_up));...
                mean(Core_2{2}.Data.Density(ind_mid));...
                mean(Core_2{2}.Data.Density(ind_down));...
                mean(Core_2{2}.Data.Density(ind_nonan))];

            OverlapPlot(Core_2, [1 2], ...
                'PlotStrat','no',...
                'lag',0,...
                'PlotDifference','no',...
                'YLimit', ylim_core,...
                'span',floor(median(thickness_act(:,ind_time) )*100));
            ylim([0 25])
            xlim([300 800])
            if count == 1
                %             yticksave = get(gca,'YTickLabel');
                %             set(gca,'YTickLabel',[])
                %             set(gca,'YTickLabel',yticksave);

                xlabel('Density (kg/m^3)')
                ylabel('Depth (m)')
            else
                set(gca,'YTickLabel',' ')
            end
        end

        for i = (count+1):num_plot
            set(f,'CurrentAxes',ha(i))
            set(gca,'Visible','off')
        end

        i_file = 1;
                NameFile = sprintf('%s/CorePlot_%s_%i.tif',OutputFolder,station{ii}, i_file);
        while exist(NameFile, 'file') == 2
            i_file = i_file + 1;
                NameFile = sprintf('%s/CorePlot_%s_%i.tif',OutputFolder,station{ii}, i_file);
        end
        print(f,NameFile,'-dpdf');

        if strcmp(vis,'off')
            close(f);
        end
    end
    ind = sum(density_meas{ii} == 0,1)==4;
    density_meas{ii}(:,ind) = [];
    density_mod{ii}(:,ind) = [];
end


end

dens_meas_all = [];
dens_mod_all = [];

for ii =1:9
    dens_meas_all = [dens_meas_all, density_meas{ii}];
    dens_mod_all = [dens_mod_all, density_mod{ii}];
end
RMSE = sqrt(nanmean((dens_meas_all-dens_mod_all).^2,2));
ME = nanmean(dens_mod_all-dens_meas_all,2);
disp('RMSE:')
disp([RMSE sum(~isnan(dens_meas_all),2)])

disp('ME:')
disp([ME sum(~isnan(dens_meas_all),2)])

lm=cell(1,9);
disp('R2:')
for i=1:4
    lm{i}=fitlm(dens_meas_all(i,:),dens_mod_all(i,:));
    disp(lm{i}.Rsquared.Adjusted);
end

%% Fig. 4: rho20 and validation with cores
f=figure('Visible',vis,'outerposition',[1 -1  30 15],'Renderer','Painter','PaperSize',[32 16]);
%  ha = tight_subplot(3,3,[0.02 0.05],[0.15 0.02],[0.07 0.1]);
ha = tight_subplot(3,2,0.07,[0.2 0.1] ,[0.1 0.01]);
col = linspecer(9);
step =48;
station3{2}= 'Dye-2';
station3{6}= 'South Dome';
station3{9}= 'Tunu-N';

for k = 1:3
set(f,'CurrentAxes',ha((k-1)*2 +1 ))
    for ii = (k-1)*3 +1:(k-1)*3 +3
        hold on
        ind = find(time_mod{ii}==datenum(2003,6,1));
        plot(time_mod{ii}(1:step:end), rho_avg_20{ii}(1:step:end)/rho_avg_20{ii}(ind)*100,'LineWidth',2,'Color',col(ii,:))
    end
    axis tight
    ylim([95 110])
    ylabel('Top 20 m average firn density \newline(% relative to 1 June 2003)','Interpreter','tex');
    box on
    legend({station3{(k-1)*3 +1:(k-1)*3 +3}},'Location','EastOutside')
%     set_monthly_tick(datenum(1998:0.001:2018,1,1))
        set(gca,'XTick',datenum(2000:5:2020,1,1),...
            'XMinorTick','on','TickLength',[0.015 0.06],...
            'XGrid','on','YGrid','on')
    datetick('x','keepticks','keeplimits')
    
    legend boxoff
    if k ~= 3
        set(gca,'XTickLabel',"")
    end
    if k ~= 2
        ylabel("")
    end
    switch k
        case 1
            text(datenum(2000,1,1),108, 'a','FontSize',15,'FontWeight','bold')
        case 2
            text(datenum(1999,1,1),108, 'b','FontSize',15,'FontWeight','bold')
        case 3
            text(datenum(1999,1,1),108, 'c','FontSize',15,'FontWeight','bold')
    end
end
ha(4).Visible = 'off';
ha(6).Visible = 'off';
ha(1).Position(3:4) = [0.27 0.22];
ha(3).Position(3:4) = [0.27 0.22];
ha(5).Position(3:4) = [0.27 0.22];
ha(2).Position(2) =0.1;
ha(2).Position(4) =1;
set(f,'CurrentAxes',ha(2))
    %  density comp overview
    l = dir(sprintf('%s/*_density_comp.txt',OutputFolder));
    path={l.folder};
    filename={l.name};
    symbol = {'<','>','v','^','d','h','x','s','x'};

    hold on
    h=[];
    diff  =[];
    plot([300 800], [300 800], 'k')
    MarkerSize = 10;
    for ii=1:length(station)
        disp(sum(~isnan(density_mod{ii}(1,:))))
        fprintf('%s\t0-5m\t%0.2f\t%0.2f\n',...
            (station{ii}),...
            (nanmean(density_meas{ii}(1,:)-density_mod{ii}(1,:))),...
            (sqrt(nanmean((density_meas{ii}(1,:)-density_mod{ii}(1,:)).^2))));
        fprintf('%s\t5-20m\t%0.2f\t%0.2f\n',...
            (station{ii}),...
            (nanmean(density_meas{ii}(2,:)-density_mod{ii}(2,:))),...
            (sqrt(nanmean((density_meas{ii}(2,:)-density_mod{ii}(2,:)).^2))));
        fprintf('%s\t>20m\t%0.2f\t%0.2f\n',...
            (station{ii}),...
            (nanmean(density_meas{ii}(3,:)-density_mod{ii}(3,:))),...
            (sqrt(nanmean((density_meas{ii}(3,:)-density_mod{ii}(3,:)).^2))));
        
        
        plot(density_meas{ii}(1,:),density_mod{ii}(1,:),...
            '.','Marker',symbol{ii},'Color',RGB('light light gray'),...
            'MarkerFaceColor',RGB('dark green'),'MarkerSize',MarkerSize);
        plot(density_meas{ii}(2,:),density_mod{ii}(2,:),...
            '.','Marker',symbol{ii},'Color',RGB('light light gray'),...
            'MarkerFaceColor',RGB('purple'),'MarkerSize',MarkerSize);    
        plot(density_meas{ii}(3,:),density_mod{ii}(3,:),...
            '.','Marker',symbol{ii},'Color',RGB('light light gray'),...
            'MarkerFaceColor',RGB('orange'),'MarkerSize',MarkerSize);
        h(ii) = plot(density_meas{ii}(1,:),NaN*density_mod{ii}(1,:),...
            '.','Marker',symbol{ii},'Color',RGB('light light gray'),...
            'MarkerFaceColor',RGB('black'),'MarkerSize',MarkerSize);
    end

    h(ii+1) = plot(density_meas{ii}(1,:),NaN*density_mod{ii}(1,:),...
        '.','Marker','o','Color',RGB('light light gray'),...
        'MarkerFaceColor',RGB('dark green'),'MarkerSize',MarkerSize);
    h(ii+2) = plot(density_meas{ii}(1,:),NaN*density_mod{ii}(1,:),...
        '.','Marker','o','Color',RGB('light light gray'),...
        'MarkerFaceColor',RGB('purple'),'MarkerSize',MarkerSize);
    h(ii+3) = plot(density_meas{ii}(1,:),NaN*density_mod{ii}(1,:),...
        '.','Marker','o','Color',RGB('light light gray'),...
        'MarkerFaceColor',RGB('orange'),'MarkerSize',MarkerSize);

    axis tight square
    xlimit=get(gca,'XLim');
     h2 = patch(xlimit([1 2 2 1]), ...
            xlimit([1 2 2 1])+[40 40 -40 -40],...
            RGB('light light gray'));
    h2.LineStyle = 'none';
    uistack(h2,'bottom')

    box on
    set(gca,'XMinorTick','on','YMinorTick','on','layer','top','YLim',xlimit)
    xlabel(sprintf('Average observed density (kg m^{-3})'),'Interpreter','tex')
    ylabel(sprintf('Average modelled density (kg m^{-3})'),'Interpreter','tex')
    legend(h([1:6 8 (end-2):end]),'CP1','Dye-2','NASA-SE','NASA-U',...
        'Saddle','South Dome','Summit',...
        '0 - 5 m depth', '5 - 20 m depth','> 20 m depth',...
        'Location','EastOutside')
    text(330,770, 'd','FontSize',15,'FontWeight','bold')

    orient('landscape')
    print(f, sprintf('%s/JOG-19-0112.Figure4.pdf',OutputFolder), '-dpdf','-bestfit')

%% Loading SEB & SMB
disp('Loading SEB files')
tic

SEB_hour = cell(1,9);
SEB_year = cell(1,9);
SEB_JJA = cell(1,9);
snow_melt = cell(1,9);
snow_melt_yearly = cell(1,9);
for i = 1:length(path_list)
   fprintf('%i/%i\n',i,length(path_list));
    % extract run parameters
    load(strcat(path_list{i},'/run_param.mat'))
    c.OutputFolder = path_list{i};
    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',path_list{i},1);

    try finfo = ncinfo(namefile);
    catch me
        namefile = sprintf('%s/%s_surf-bin-%i.nc',path_list{i},station{i},1);
        finfo = ncinfo(namefile);
    end
    
    names={finfo.Variables.Name};
    for ii= 1:size(finfo.Variables,2)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{ii}), namefile,char(names{ii})));
    end
    time_mod2 = datenum(Year,1,Day,Hour,0,0);
    if ~isfield(c,'verbose')
        c.verbose = 1;
    end
    % extracting observed
    SEB_hour{i} = table(time_mod2,SHF.*c.dt_obs,LHF,SRin.*c.dt_obs , SRout.*c.dt_obs,LRin.*c.dt_obs, LRout_mdl.*c.dt_obs,...
        rainHF.*c.dt_obs,GF.*c.dt_obs, meltflux.*c.dt_obs, meltflux.*c.dt_obs./c.dev./c.L_fus./c.rho_water,...
        'VariableNames',{'time','SHF_J','LHF_J','SRin_J','SRout_J','LRin_J','LRout_J','RainHF_J','GF_J','MeltEnergy_J',...
        'Melt_mweq'});
    SEB_hour{i}.SRnet_J = SEB_hour{i}.SRin_J-SEB_hour{i}.SRout_J;
    SEB_hour{i}.LRnet_J = SEB_hour{i}.LRin_J-SEB_hour{i}.LRout_J;
    SEB_JJA{i} = AvgTableJJA(SEB_hour{i},'sum');
    SEB_year{i} = AvgTable(SEB_hour{i},'yearly','sum');

% Loading SMB
    snow_melt{i} = table(time_mod2, snowfall, meltflux/334*3600/1000000, sublimation, ...
        'VariableNames',{'time','Snowfall_mweq','Melt_mweq','Sublimation_mweq'});
    snow_melt_yearly{i} = AvgTable(snow_melt{i},'yearly','sum');
end

%% Validation accumulation
time_accum =1990:2017;
%  loading accumulation records
opts = delimitedTextImportOptions("NumVariables", 5);
opts.DataLines = [2, Inf];
opts.Delimiter = ";";
opts.VariableNames = ["ID", "Name", "Latitude", "Longitude", "Citation"];
opts.VariableTypes = ["double", "string", "double", "double", "categorical"];
opts = setvaropts(opts, 2, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 5], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
metadata = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Code\Accumulation\metadata.csv", opts);
clear opts
s = char(97:122);

f = figure('Visible',vis,'Outerposition',[1 1 23 18]);
ha = tight_subplot(3,3,[0.03 0.015],[0.15 0.01],[0.07 0.22]);
col = linspecer(15,'sequential');
% col = col(randperm(size(col,1)),:);
count = 1;
    leg_text{1} = 'Station-derived';
clearvars h
for ii =1:length(station) 
    load(strcat(path_list{ii},'/run_param.mat'))

    disp(' ')
    for i = 1:size(metadata,1)

        switch metadata.Name{i}
            case {'PARCA_6642','PARCA_nasau','PARCA_saddlea','PARCA_nasaea','Site 10'}
                metadata.Latitude(i)=0;
            case 'PARCA-1998-cores (CORE 6642 (B))'
                metadata.Name{i}='Core 6642-B';

            case 'Summit-Zoe-10'
                metadata.Name{i}='Owen';
        end
        if ~isempty(strfind( metadata.Name{i},'PARCA'))
    %         disp(metadata.Name{i})
                metadata.Latitude(i)=0;
        end
    end

    disp(c.station)
    ind = distance(c.lat,c.lon,metadata.Latitude,metadata.Longitude)<0.05;
    if ii == 6
            ind = distance(c.lat,c.lon,metadata.Latitude,metadata.Longitude)<0.3;
    end

    ind=find(ind);
    year=cell(1,9);
    accum={};
    counter = 1;
    ind_remove = [];
    for i = 1:length(ind)
        opts = delimitedTextImportOptions("NumVariables", 2);
        opts.DataLines = [1, Inf];
        opts.Delimiter = ";";
        opts.VariableNames = ["year", "b_mm"];
        opts.VariableTypes = ["double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        accum{counter} = readtable(sprintf('../Accumulation/data/%i.csv',ind(i)), opts);
        clear opts
        
        if max(accum{counter}.year)<1990
            ind_remove = [ind_remove ind(i)];
            continue
        end        
        counter=counter+1;
    end
        accum(ismember(ind,ind_remove)) = [];
        ind(ismember(ind,ind_remove)) = [];

    DV = datevec(snow_melt_yearly{ii}.time);
    years = DV(:,1);
    SMB_station = snow_melt_yearly{ii}.Snowfall_mweq ...
        -snow_melt_yearly{ii}.Sublimation_mweq;
    
    set(f,'CurrentAxes',ha(ii))
    hold on   
    h(1) = plot(years,SMB_station,':ok','LineWidth',2,...
        'MarkerFaceColor','k','MarkerSize',4);
    lm = fitlm(years,SMB_station);
    plot(years,...
        lm.Coefficients.Estimate(1) + lm.Coefficients.Estimate(2) .* years,...
        'm','LineWidth',1.5)
    if ii == 9
        text(1993,0.3,sprintf('trend: %+0.0f mm weq dec^{-1}',...
        lm.Coefficients.Estimate(2)*10*1000),...
        'Interpreter','tex',...
        'Color','m',...
        'FontWeight','bold','FontSize',10);
    else
        text(1993,0.18*max(get(gca,'YLim')),sprintf('trend: %+0.0f mm weq dec^{-1}',...
        lm.Coefficients.Estimate(2)*10*1000),...
        'Interpreter','tex',...
        'Color','m',...
        'FontWeight','bold','FontSize',10);
    end
    
    for i =1:length(accum)
        leg_text{count+1} = metadata.Name{ind(i)};
        [year_ordered, ind_o] = sort(accum{i}.year);
        h(count+1) = plot(year_ordered, accum{i}.b_mm(ind_o),...
            ':o','LineWidth',2,'Color',col(count,:),...
            'MarkerFaceColor',col(count,:),'MarkerSize',4);
        count = count+1;
    end  
 
    uistack(h(1),'top');
    xlim(time_accum([1 end]))
    set(gca,'layer','top')

    switch ii
        case {1,2,3}
           set(gca,'YLim',[0 1.3])
           text(time_accum(1)+1,1.2,[s(ii) ' ' station{ii}],'FontSize',15)
        case {4,5,6}
           set(gca,'YLim',[0 1])
           text(time_accum(1)+1,0.9,[s(ii) ' ' station{ii}],'FontSize',15)
        case {7,8,9}
            set(gca,'YLim',[0 .5])
           text(time_accum(1)+1,0.45,[s(ii) ' ' station{ii}],'FontSize',15)
    end
   set (gca,'YMinorTick','on','XMinorTick','on',...
       'XTick',time_accum(1:5:end),...
       'TickLength', [0.01 0.025].*2)
    if ismember(ii,1:6)
       set(gca,'XTickLabel','')
    else
        set(gca,'XTickLabel',time_accum(1:5:end)','XTickLabelRotation',45)
   end
   if ismember(ii,[2 3 5 6 8 9])
       set(gca,'YTickLabel','')
   end
   if ii == 4
       h_label = ylabel( sprintf('SMB (m w.e.)'),'Interpreter','tex');
   end
   if ii == 8
       xlabel('Year')
   end
   
   box on

end

legendflex(h,leg_text, 'ref', gcf, ...
                       'anchor', {'e','e'}, ...
                       'buffer',[-10 30], ...
                       'ncol',1, ...
                       'fontsize',11,...
                       'box','off',...
                       'Interpreter','none');
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
ppos(3:4) = pos(3:4);
% pos(1:2) = [1 1];
set(gcf,'paperposition',ppos);
set(gcf,'units',unis);
print(f, sprintf('%s/Accumulation_val',OutputFolder), '-dpdf','-bestfit')

%% Fig. 3: Validation accumulation (gray version)

time_accum =1990:2017;

%  loading accumulation records
opts = delimitedTextImportOptions("NumVariables", 5);
opts.DataLines = [2, Inf];
opts.Delimiter = ";";
opts.VariableNames = ["ID", "Name", "Latitude", "Longitude", "Citation"];
opts.VariableTypes = ["double", "string", "double", "double", "categorical"];
opts = setvaropts(opts, 2, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 5], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
metadata = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Code\Accumulation\metadata.csv", opts);
clear opts
s = char(97:122);

f = figure('Visible',vis,'Outerposition',[1 1 20 18],'Renderer','Painter','PaperSize',[22 16]);
ha = tight_subplot(3,3,0.012,[0.15 0.1],[0.11 0.02]);

count = 1;
    leg_text{1} = 'Station-derived';
clearvars h
for ii =1:length(station) 
    load(strcat(path_list{ii},'/run_param.mat'))

    pit_data = [];
    if size(pit_data,1)>0  
        leg_text{2} = 'Snowpit-derived';
        count = count+1;

        pit_data.SMB_sta = NaN(size(pit_data,1),1);
        pit_data.tot_subl = NaN(size(pit_data,1),1);

        for i=1:size(pit_data,1)
            % closest time step to the snowpit survey date
            date_end = datenum(pit_data.Date(i,:));
            [~, pit_data.ind_end(i)] = min(abs(date_end - SEB_hour{ii}.time));
            [~, pit_data.ind_end_in_subl(i)] = min(abs(date_end - subl.time));
            temp = datevec(date_end);
        
            % We assume winter accumulation starts on 1st September
            date_start = datenum(temp(:,1)-1,09,01);
            [~, pit_data.ind_start(i)] = min(abs(SEB_hour{ii}.time-date_start));
            [~, pit_data.ind_start_in_subl(i)] = min(abs(date_start - subl.time));
            date_start = SEB_hour{ii}.time(pit_data.ind_start(i));
        
            % summing the sublimation that occured since the 1st Sept. until
            % the snowpit survey
            pit_data.tot_subl(i) = sum(subl.estim( ...
                pit_data.ind_start_in_subl(i):pit_data.ind_end_in_subl(i))); % in mm weq

            pit_data.SMB_sta(i) = sum(snow_melt{ii}.Snowfall_mweq(...
                pit_data.ind_start(i):pit_data.ind_end(i))) ...
                + pit_data.tot_subl(i);
            if pit_data.ind_end(i) ~= pit_data.ind_start(i)
                missed_sf = sum(snow_melt{ii}.Snowfall_mweq(...
                    pit_data.ind_end(i):...
                    min(pit_data.ind_start(i)+365*24,length(snow_melt{ii}.Snowfall_mweq))))*1000;
                pit_data.SWE_pit_cor(i) = pit_data.SWE_pit(i) + missed_sf;
                fprintf('%0.1f mm weq snowfall was monitored after the snow pit\n',missed_sf);
            else
                pit_data.SWE_pit_cor(i) = pit_data.SWE_pit(i);
            end
            
        end
%             pit_data(pit_data.SMB_sta<0.1,:)=[];

    end

    disp(' ')

for i = 1:size(metadata,1)

    switch metadata.Name{i}
        case {'PARCA_6642','PARCA_nasau','PARCA_saddlea','PARCA_nasaea','Site 10'}
            metadata.Latitude(i)=0;
        case 'PARCA-1998-cores (CORE 6642 (B))'
            metadata.Name{i}='Core 6642-B';

        case 'Summit-Zoe-10'
            metadata.Name{i}='Owen';
    end
    if ~isempty(strfind( metadata.Name{i},'PARCA'))
%         disp(metadata.Name{i})
            metadata.Latitude(i)=0;
    end
end
    disp(c.station)
    ind = distance(c.lat,c.lon,metadata.Latitude,metadata.Longitude)<0.05;
    if ii == 6
        ind = distance(c.lat,c.lon,metadata.Latitude,metadata.Longitude)<0.3;
    end

    ind=find(ind);
    year=cell(1,9);
    accum={};
    counter = 1;
ind_remove = [];
    for i = 1:length(ind)
        opts = delimitedTextImportOptions("NumVariables", 2);

        opts.DataLines = [1, Inf];
        opts.Delimiter = ";";

        opts.VariableNames = ["year", "b_mm"];
        opts.VariableTypes = ["double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        accum{counter} = readtable(sprintf('../Accumulation/data/%i.csv',ind(i)), opts);
        clear opts
        
        if max(accum{counter}.year)<1990
            ind_remove = [ind_remove ind(i)];
            continue
        end        
        counter=counter+1;
    end
        accum(ismember(ind,ind_remove)) = [];
        ind(ismember(ind,ind_remove)) = [];

    DV = datevec(snow_melt_yearly{ii}.time);
    years = DV(:,1);
    SMB_station = snow_melt_yearly{ii}.Snowfall_mweq ...
        -snow_melt_yearly{ii}.Sublimation_mweq;
    
    set(f,'CurrentAxes',ha(ii))
    hold on   
    h(1) = plot(years,SMB_station*1000,':ok','LineWidth',2,...
        'MarkerFaceColor','k','MarkerSize',4);
    lm = fitlm(years,SMB_station*1000);
    h(3) = plot(years,...
        lm.Coefficients.Estimate(1) + lm.Coefficients.Estimate(2) .* years,...
        'Color',RGB('red'),'LineWidth',2.5);
    disp(max(lm.Coefficients.pValue))
   switch ii
      case {7,8,9}
        ylim([0 500])
       otherwise
        ylim([0 1400])
   end
    if ismember(ii , [3 7 8])
        text(1991,0.1*max(get(gca,'YLim')),sprintf('trend: %+0.0f mm w.e. decade^{-1}',...
        lm.Coefficients.Estimate(2)*10),...
        'Interpreter','tex',...
        'Color',RGB('red'),...
        'FontWeight','bold','FontSize',10.5);
    else
        text(1991,0.7*max(get(gca,'YLim')),sprintf('trend: %+0.0f mm w.e. decade^{-1}',...
        lm.Coefficients.Estimate(2)*10),...
        'Interpreter','tex',...
        'Color',RGB('red'),...
        'FontWeight','bold','FontSize',10.5);
    end
    
    for i =1:length(accum)
        [year_ordered, ind_o] = sort(accum{i}.year);
        h(count+1) = plot(year_ordered, accum{i}.b_mm(ind_o)*1000,...
            ':o','LineWidth',2,'Color',RGB('dark gray'),...
            'MarkerFaceColor',RGB('dark gray'),'MarkerSize',4);
        count = count+1;
    end  
 leg_text = {'Derived from the station',...
     'Derived from firn and ice cores',...
     'Linear trend in station-derived accumulation'};
    uistack(h(1),'top');
    
    xlim(time_accum([1 end]))
    set(gca,'layer','top')
    switch ii
        case {2}
            tmp = [s(ii) ' Dye-2' ];
         case {6}
            tmp = [s(ii) ' South Dome' ];       
        case {9}
            tmp = [s(ii) ' Tunu-N'];
        otherwise
            tmp = [s(ii) ' ' station{ii}];
    end
    text(time_accum(1)+1,0.9*max(get(gca,'YLim')),tmp,'FontSize',15,'FontWeight','bold')

   set (gca,'YMinorTick','on','XMinorTick','on',...
       'XTick',time_accum(1:5:end),...
       'TickLength', [0.01 0.025].*2)
    if ismember(ii,1:6)
       set(gca,'XTickLabel','')
    else
        set(gca,'XTickLabel',time_accum(1:5:end)','XTickLabelRotation',45)
   end
   if ismember(ii,[2 3 5 6 8 9])
       set(gca,'YTickLabel','')
   end
   if ii == 4
       h_label = ylabel( sprintf('Snow accumulation (mm w.e.)'),'Interpreter','tex');
   end
   if ii == 8
       xlabel('Year')
   end
   set (gca,'YGrid','on','XGrid','on')
   box on
end

legendflex(h(1:3),leg_text, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 7], ...
                       'ncol',1, ...
                       'fontsize',11,...
                       'box','off',...
                       'Interpreter','none');
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
ppos(3:4) = pos(3:4);
% pos(1:2) = [1 1];
set(gcf,'paperposition',ppos);
set(gcf,'units',unis);
print(f, sprintf('%s/JOG-19-0112.Figure3.pdf',OutputFolder), '-dpdf','-bestfit')

%% Fig. 6: Firn heat budget 
f=figure('Visible',vis,'outerposition',[1 -1  25 25],'Renderer','Painter','PaperSize',[22 19]);
ha = tight_subplot(3,3,[0.02 0.02],[0.17 0.1],[0.09 0.1]);

 time_mod=cell(1,9);
for ii = 1:length(path_list)
   
     set(f,'CurrentAxes',ha(ii))
    filename = [path_list{ii} '/track_CC.nc'];
%     disp(filename)
    finfo  = ncinfo(filename);
    CC20 = ncread(filename,'CC20');
    time_mod{ii} = datenum(1998,6,1):1/24:(datenum(1998,6,1)+length(CC20(1,:))/24-1/24);
    
    dCC_snow = CC20(2,:)-CC20(1,:);
    dCC_subl = CC20(3,:)-CC20(2,:);
    dCC_melt = CC20(4,:)-CC20(3,:);
    dCC_rfrz = CC20(6,:)-CC20(5,:);
    dCC_SEB = CC20(6,:)*0;
    dCC_SEB(2:end) = CC20(1,2:end)- CC20(6,1:end-1);
    
     dCC = table(time_mod{ii}',...
         dCC_subl', dCC_melt', dCC_snow',dCC_rfrz',dCC_SEB',...
         'VariableNames',{'time','subl','melt','snow','rfrz','SEB'});
    dCC_avg = AvgTable(dCC,'daily','sum');
    
    CC_track = table(time_mod{ii}',CC20(end,:)',...
         'VariableNames',{'time','CC20'});
    CC_avg = AvgTable(CC_track,'daily','mean');

    DV  = datevec(dCC_avg.time);  % [N x 6] array
    DV  = DV(:, 1:3);   % [N x 3] array, no time
    DV2 = DV;
    DV2(:, 2:3) = 0;    % [N x 3], day before 01.Jan
    DOY = datenum(DV) - datenum(DV2);
    
    CC_climatology = array2table([1:366]','VariableNames',{'doy'});
    for i=1:366
        CC_climatology.dCC_snow(i) = mean(dCC_avg.snow(DOY==i));
        CC_climatology.dCC_subl(i) = mean(dCC_avg.subl(DOY==i));
        CC_climatology.dCC_melt(i) = mean(dCC_avg.melt(DOY==i));
        CC_climatology.dCC_rfrz(i) = mean(dCC_avg.rfrz(DOY==i));
        CC_climatology.dCC_SEB(i) = mean(dCC_avg.SEB(DOY==i));
        CC_climatology.CC20(i) = mean(CC_avg.CC20(DOY==i));
    end
    
    [AX,H1,H2] = plotyy(1:366,(CC_climatology.dCC_SEB+...
        CC_climatology.dCC_snow +...
        CC_climatology.dCC_rfrz)./1000,...
        1:366, CC_climatology.CC20*100./CC_climatology.CC20(1));

    H1.LineWidth=0.0000001;
    H2.Color=RGB('dark gray');
    H1.Color=RGB('light light gray');
    AX(2).YAxis.Color    =H2.Color;
    AX(1).YAxis.Color    ='k';
    H2.LineWidth=2;
    hold on
    axis(AX(1))
    hold on

    
    aux=CC_climatology.dCC_SEB;
    aux(aux>0)=0;
    hh_d=area(1:366, (aux+CC_climatology.dCC_snow+CC_climatology.dCC_rfrz)./1000,...
        'LineStyle','none');
    hh2=area(1:366,(CC_climatology.dCC_rfrz+CC_climatology.dCC_snow)./1000,...
        'LineStyle','none');
    hh3=area(1:366,(CC_climatology.dCC_snow)./1000,...
        'LineStyle','none');
    
    aux=CC_climatology.dCC_SEB;
    aux(aux<0)=0;
    hh1=area(1:366,(aux + CC_climatology.dCC_melt+CC_climatology.dCC_subl)./1000,...
        'LineStyle','none');
    col = hh1.FaceColor;
    hh1.FaceColor = hh_d.FaceColor;

    hh4 = area(1:366,(CC_climatology.dCC_melt+CC_climatology.dCC_subl)./1000,...
        'LineStyle','none');
    hh4.FaceColor = col;


    ylim([-1.8 1.2])
    AX(1).YTick=-2:1:2;
    box off
    xlim([1 365])
    AX(1).XTick=[1 61 122 183  245 306];
    AX(2).XTick=[1 61 122 183  245 306];
    AX(2).XTickLabel = {'Jan.' 'Mar.'  'May' 'Jul.' 'Sept.'  'Nov.' };
     AX(1).XMinorTick = 'on';
     AX(2).XMinorTick = 'on';
    AX(1).XAxis.MinorTickValues=[1 31 61 92 122 153 183 214 245 275 306 336];
    AX(2).XAxis.MinorTickValues=[1 31 61 92 122 153 183 214 245 275 306 336];
     AX(1).YMinorTick = 'on';
     AX(2).YMinorTick = 'on';
     
     axes(AX(2))
     ylim([92.5 105])
    AX(2).YTick=90:5:105;
    s = char(97:122);
    h_leg = legend('snosssssssssw \newline ','Interpreter','tex','Location','SouthWest');
    

      h_text =  annotation('textbox',...
            h_leg.Position+[0 0.01 0 0],...
            'String',[s(ii) ' ' station3{ii}],'FontSize',15,'FontWeight','bold','LineStyle','none');
      h_leg.Visible =  'off';

      if ismember(ii,[2 3 5 6 8 9])
          AX(1).YTickLabel = '';
      end
      if ismember(ii,[1 2 4 5 7 8])
          AX(2).YTickLabel = '';
      end
      if ismember(ii,1:6)
          AX(2).XTickLabel = '';
      end
          AX(1).XTickLabel = '';
          
      if ii ==  4
          ylabel(AX(1), 'Daily contribution to cold content (MJ m^{-2} d^{-1})','Interpreter','tex');
      elseif ii==6
         ylabel( AX(2),'Cold content (% relative to 1st Jan.)');
      elseif ii==8
          xlabel('Month')
      end
          AX(1).TickLength=AX(1).TickLength*4;
      AX(2).TickLength=AX(2).TickLength*4;
      AX(2).XTickLabelRotation=45;
      AX(1).XTickLabelRotation=45;
       AX(1).XAxisLocation = 'top';
      AX(2).XAxis.Visible = 'on';
      AX(2).XGrid='on';
      box off
    xlim([1 365])
end
legendflex([hh1, hh2, hh3, hh4, H2],...
    {'heat conduction','latent-heat release', 'accumulation','ablation','cold content'},...
    'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[-10 0], ...
                       'nrow',2, ...
                       'fontsize',15,...
                       'box','off',...
                       'Interpreter','none');
                   
print(f, sprintf('%s/JOG-19-0112.Figure6.pdf',OutputFolder), '-dpdf','-bestfit')

%% Fig. 7: Firn heat budget cumulative
f=figure('Visible',vis,'outerposition',[1 -1  25 25],'Renderer','Painter','PaperSize',[25 20]);
 ha = tight_subplot(3,3,0.03,[0.19 0.1], 0.11);

time_mod=cell(1,9);
cc_summary=cell(1,9);
cc_summary{1} = sprintf('station;SEB;SEB_y;snow;snow_y;rfrz;rfrz_y;subl;subl_y;melt;melt_y\n');

for ii = 1:length(path_list)
   
     set(f,'CurrentAxes',ha(ii))
    filename = [path_list{ii} '/track_CC.nc'];
%     disp(filename)
    finfo  = ncinfo(filename);
    CC20 = ncread(filename,'CC20');
    time_mod{ii} = datenum(1998,6,1):1/24:(datenum(1998,6,1)+length(CC20(1,:))/24-1/24);

    %             CC20(1,k) = c.rhoCC20_aft_comp(2);
    %             CC20(2,k) = c.rhoCC20_aft_snow(2);
    %             CC20(3,k) = c.rhoCC20_aft_subl(2);
    %             CC20(4,k) = c.rhoCC20_aft_melt(2);
    %             CC20(5,k) = c.rhoCC20_aft_runoff(2);
    %             CC20(6,k) = c.rhoCC20_aft_rfrz(2);
    dCC_snow = CC20(2,:)-CC20(1,:);
    dCC_subl = CC20(3,:)-CC20(2,:);
    dCC_melt = CC20(4,:)-CC20(3,:);
    dCC_rfrz = CC20(6,:)-CC20(5,:);
    dCC_SEB = CC20(6,:)*0;
    dCC_SEB(2:end) = CC20(1,2:end)- CC20(6,1:end-1);
    
     dCC = table(time_mod{ii}',...
         dCC_subl', dCC_melt', dCC_snow',dCC_rfrz',dCC_SEB',...
         'VariableNames',{'time','subl','melt','snow','rfrz','SEB'});
    dCC_avg = AvgTable(dCC,'daily','sum');
    
    CC_track = table(time_mod{ii}',CC20(end,:)',...
         'VariableNames',{'time','CC20'});
    CC_avg = AvgTable(CC_track,'daily','mean');

    [AX,H1,H2] = plotyy(dCC_avg.time,cumsum(dCC_avg.SEB)./1000,...
        CC_avg.time,CC_avg.CC20*100./CC_avg.CC20(...
        CC_avg.time==datenum(1998,6,1)));
%     H1 = plot(dCC_avg.time,cumsum(dCC_avg.SEB)./1000);

    H1.LineWidth=0.0000001;
    H2.Color=RGB('dark gray');
    H2.LineWidth=1;

    H1.Color=RGB('light light gray');
    AX(2).YAxis.Color    =H2.Color;
    AX(1).YAxis.Color    ='k';
    axis(AX(1))
    hold on
    hh1=plot(dCC_avg.time,cumsum(dCC_avg.SEB)./1000,...
        'LineWidth',2.5);

    hh0=plot(dCC_avg.time,cumsum(dCC_avg.SEB)*0,'--k');

    hh2=plot(dCC_avg.time,cumsum(dCC_avg.rfrz)./1000,...
        'LineWidth',2.5);
    hh3=plot(dCC_avg.time,cumsum(dCC_avg.snow)./1000,...
        'LineWidth',2.5);
    hh4=plot(dCC_avg.time,cumsum(dCC_avg.subl+dCC_avg.melt)./1000,...
        'LineWidth',2.5);

    cc_summary{1+ii} = sprintf('%s;%0.2f;%0.2f;%0.2f;%0.2f;%0.2f;%0.2f;%0.2f;%0.2f;%0.2f;%0.2f\n',...
    station{ii},...
    max(cumsum(dCC_avg.SEB))./1000,...
    max(cumsum(dCC_avg.SEB))./1000/20,...
    min(cumsum(dCC_avg.snow))./1000,...
    min(cumsum(dCC_avg.snow))./1000/20,...
    min(cumsum(dCC_avg.rfrz))./1000,...
    min(cumsum(dCC_avg.rfrz))./1000/20,...
    max(cumsum(dCC_avg.subl))./1000,...
    max(cumsum(dCC_avg.subl))./1000/20,...
    max(cumsum(dCC_avg.melt))./1000,...
    max(cumsum(dCC_avg.melt))./1000/20);    

     switch ii
      case {1,2,3}
        ylim([-2000 2000])
        AX(1).YTick=-2000:1000:2000;
      case {4,5,6}
        ylim([-1200 1200])
        AX(1).YTick=-1200:600:1200;
      case {7,8,9}
        ylim([-400 400])
        AX(1).YTick=-400:200:400;
     end
           box off

     axes(AX(2))
    ylim([70 130])
    AX(2).YTick=70:15:130;
    AX(2).YGrid='on';
    s = char(97:122);
    h_leg = legend(sprintf('snosssssssss\nw'),'Location','SouthWest');
    

      h_text =  annotation('textbox',...
            h_leg.Position+[0 0 0 0],...
            'String',[s(ii) ' ' station3{ii}],'FontSize',15,'FontWeight','bold','LineStyle','none');
      h_leg.Visible =  'off';
      
      AX(2).XTick = datenum(1990:2:2017,1,1);
      AX(1).XTick = datenum(1990:2:2017,1,1);
      AX(1).XAxis.MinorTickValues = datenum(1998:1:2017,1,1);
      AX(2).XGrid = 'on';
      datetick('x','yyyy','keepticks','keeplimits')
      xticklabels = get(gca,'XTickLabel');
        for i =1:length(xticklabels)
            if floor(i/2)==i/2
                xticklabels(i,:) = '    ';
            end
        end
       AX(2).XTickLabelRotation = 45;
      if ismember(ii,[2 3 5 6 8 9])
          AX(1).YTickLabel = '';
      end
      if ismember(ii,[1 2 4 5 7 8])
          AX(2).YTickLabel = '';
      end
      if ismember(ii,1:6)
          AX(2).XTickLabel = '';
      end
          AX(1).XTickLabel = '';
          
      if ii ==  4
          ylabel(AX(1), 'Cumulated contribution to cold content (MJ m^{-2})','Interpreter','tex');
      elseif ii==6
         ylabel( AX(2),'Cold content (% relative to June 1998)');
      elseif ii==8
          xlabel(AX(2), 'Year')
      end
      
      AX(1).TickLength=AX(1).TickLength*4;
      AX(2).TickLength=AX(2).TickLength*4;

       AX(1).XAxisLocation = 'top';
      AX(2).XAxis.Visible = 'on';
      AX(1).XAxis.MinorTickValues = datenum(1998:1:2017,1,1);
      AX(2).XAxis.MinorTickValues = datenum(1998:1:2017,1,1);
      AX(1).XMinorTick = 'on';
      AX(2).XMinorTick = 'on';

      box off
      uistack(H2,'bottom')
      
end
legendflex([hh1, hh2, hh3,hh4, H2],...
    {'heat conduction','latent-heat release', 'accumulation',...
    'ablation' ,'cold content'},...
    'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[-10 0], ...
                       'nrow',2, ...
                       'fontsize',15,...
                       'box','off',...
                       'Interpreter','none');
cc_summary = cc_summary';
namesummary=sprintf('%s/cc_summary.csv',OutputFolder);
delete(namesummary);
fid=fopen(namesummary,'a+');
for k = 1:length(cc_summary)
    fprintf(fid,cc_summary{k});
end
fclose(fid);
    
print(f, sprintf('%s/JOG-19-0112.Figure7.pdf',OutputFolder), '-dpdf','-bestfit')
% clearvars dCC dCC_avg dCC_melt dCC_month dCC_rfrz dCC_SEB dCC_snow dCC_subl cc_summary ...
%     hh1 hh2 hh3 hh4 hh5 H1 H2 AX CC20 CC_avg a aux aux_sum count fid namesummary namefile ...
%     h_leg h_text hh6 hh7 hh_d i i_core ind ind_bin ind_bot ind_CP ind_Dye2 ind_lim ...
%     ind_nan ind_nonan ind_remove ind_start ind_Summit ind_year index indexToDupes k_eff ...
%     list lm rho_20 rho_alt step s Summit_Miller depth_extra_Dye2 depth_extra_Summit ...
%     T_20 T_diff T_ice_alt temp th thickness_alt tmp tmp_text TT_CP TT_sum VarName x y xlimit ...
%     ylimit year_uni AX_d CC_track CP depth_alt depth_alt_d2 ...
%     depth_temp finfo ha h_tit ii j k kk Core

%% Loading tilt-corrected and uncorrected SEB
% error('stopping here before clearing variables')
station_list=station;
for kk =2:-1:1
    if kk==1
        WorkingFolder = './Output/Corrected';
    else
        WorkingFolder = './Output/no cor';
    end
    list = dir(WorkingFolder);
    folder_list = {list.name};
    folder_list = folder_list(3:end)';
    ind_delete = zeros(length(folder_list),1);
    for i = 1:length(folder_list)
        ind_delete(i) = 1;
        for j = 1:length(station_list)
            if (~isempty(strfind(folder_list{i},station_list{j})))&&...
                    (isempty(strfind(folder_list{i},'Erro')))
                ind_delete(i) = 0;
            end
        end
    end
    folder_list(ind_delete==1)= [];

    path_list{kk}=cell(1,9);
    station=cell(1,9);
    count = 0;
    for i = [1 2 4 5 6 7 3 8 9] % 1:length(folder_list)
        count = count +1;
        path_list{kk}{count} = [WorkingFolder '/' folder_list{i}];
        temp = folder_list{i};
        ind =strfind(temp,'_0');
        station{count} = temp(1:ind-1);
    end
end
  
station2 = station;
for i =1:length(station)
    station2{i} = strrep(station{i},'-','');
    station2{i} = strrep(station2{i},' ','');
end

% Loading SEB a files if needed
disp('Loading SEB files')
tic
origin_table_1 = table;
origin_table_2 = table;
for kk =1:2
    SEB_hour{kk} = cell(1,9);
    SEB_year{kk} = cell(1,9);
    SEB_JJA{kk} = cell(1,9);
    for i = 1:length(path_list{kk})

       fprintf('%i/%i\n',i,length(path_list{kk}));

        file_seb_year = sprintf('%s/SEB_year.txt',path_list{kk}{i});
        file_seb_JJA = sprintf('%s/SEB_JJA.txt',path_list{kk}{i});
    
%     if or(~exist(file_seb_year,'file'),~exist(file_seb_JJA,'file'))
        % extract run parameters
        load(strcat(path_list{kk}{i},'/run_param.mat'))
        c.OutputFolder = path_list{kk}{i};
        % extract surface variables
        namefile = sprintf('%s/surf-bin-%i.nc',path_list{kk}{i},1);
        try finfo = ncinfo(namefile);
        catch me
            namefile = sprintf('%s/%s_surf-bin-%i.nc',path_list{kk}{i},station{i},1);
            finfo = ncinfo(namefile);
        end
        names={finfo.Variables.Name};
        for ii= 1:size(finfo.Variables,2)
            eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{ii}), namefile,char(names{ii})));
        end
        time_mod = datenum(Year,1,Day,Hour,0,0);
        if ~isfield(c,'verbose')
            c.verbose = 1;
        end
        % extracting observed
    c.InputAWSFile = ['./Output/Corrected/AWS data/' c.InputAWSFile(26:end)];
            [~, ~, ~, ~, ~,...
        ~, ~, ~, ~, ~,~, ...
        ~, ~, ~, ~, ~, ~, ...
        ~, ~, ~, ~, ~, ~,...
        ~, ~, LRin, ~, ~, ...
        ~, ~, ~, ~, c] = ...
        ExtractAWSData(c);

    if exist('SWin','var')
        SRin = SWin;
    end
    if exist('SWout','var')
        SRout = SWout;
    end
        SEB_hour{kk}{i} = table(time_mod,SHF,LHF,SRin , SRout,LRin, LRout_mdl,...
            rainHF,GF, meltflux, meltflux.*c.dt_obs./c.dev./c.L_fus./c.rho_water,...
            'VariableNames',{'time','SHF_Wm2','LHF_Wm2','SRin_Wm2','SRout_Wm2','LRin_Wm2','LRout_Wm2','RainHF_Wm2','GF_Wm2','MeltEnergy_Wm2',...
            'Melt_mweq'});
        SEB_hour{kk}{i}.SRnet_Wm2 = SEB_hour{kk}{i}.SRin_Wm2-SEB_hour{kk}{i}.SRout_Wm2;
        SEB_hour{kk}{i}.LRnet_Wm2 = SEB_hour{kk}{i}.LRin_Wm2-SEB_hour{kk}{i}.LRout_Wm2;
        SEB_JJA{kk}{i} = AvgTableJJA(SEB_hour{kk}{i},'mean');
        SEB_year{kk}{i} = AvgTable(SEB_hour{kk}{i},'yearly2','mean');
        tmp1 = AvgTableJJA(SEB_hour{kk}{i},'sum');
        SEB_JJA{kk}{i}.Melt_mweq  = tmp1.Melt_mweq;
        tmp2 = AvgTable(SEB_hour{kk}{i},'yearly2','sum');
        SEB_year{kk}{i}.Melt_mweq  = tmp2.Melt_mweq;
    end

    for ii = 1:length(names)
        clearvars(names{ii})
    end
end
    
% writetableFast(origin_table_1, sprintf('%s/origin_table_1.csv',OutputFolder),';') 
toc

%% Fig. 2: Tilt-corrected melt
f = figure('Visible',vis,'outerposition',[0 0 25 18],'Renderer','Painter','PaperSize',[25 19]);
ha = tight_subplot(3,3,0.043,[0.1 0.08],0.07);
for ii =1:length(path_list{1})
    set(f,'CurrentAxes',ha(ii))
    hold on
    melt_cor = SEB_JJA{1}{ii}.Melt_mweq;
    melt_ucor = SEB_JJA{2}{ii}.Melt_mweq;
    time_mod = SEB_JJA{1}{ii}.time;
    tt_cor = sum(melt_cor)*1000;
    tt_ucor = sum(melt_ucor)*1000;
    
    stairs([time_mod;  time_mod(end)+365], ...
        [melt_cor; melt_cor(end)].*1000,...
        'LineWidth',2.5)
    stairs([time_mod;  time_mod(end)+365], ...
        [melt_ucor; melt_ucor(end)].*1000,...
       'LineWidth',2)
%     plot([time_mod(1)  time_mod(end)+365] ,[0 0],'--k')
t0 = text(datenum(1999,1,1), 1.1*max(melt_cor)*1000*0.81,'Total:');
t0.FontWeight = 'bold';
col = lines(2);
t1 = text(datenum(1999,1,1), 1.1*max(melt_cor)*1000*0.7,sprintf('%0.1f mm',tt_cor));
t1.FontWeight = 'bold';
t1.Color = col(1,:);
t2 = text(datenum(1999,1,1), 1.1*max(melt_cor)*1000*0.59,sprintf('%0.1f mm',tt_ucor));
t2.Color = col(2,:);
t2.FontWeight = 'bold';

    h_text = set_plot_style(ii, station, time_mod, ...
        'Annual melt (mm)');
    ylimits = get(gca,'YLim');
    set(gca,'YLim',[ylimits(1) ylimits(2)*1.2]);
    h_text.Units = 'Normalized';
    h_text.Position(1:2) = [0.03 0.9];
    ha(ii).XMinorTick='on';
    ha(ii).XAxis.MinorTickValuesMode='auto';
    ha(ii).TickLength = ha(ii).TickLength*2;
    if ii == 7
        t0.Position(2) =  35*0.75;
        t1.Position(2) =  35*0.64;
        t2.Position(2) =  35*0.53;

        set(gca,'YLim',[ylimits(1) 35]);
    end
end
legendflex({'applying tilt correction on radiation data', ...
    'without tilt correction on the radiation data'},...
    'ref',gcf,...
    'anchor',{'n','n'},...
    'buffer',[0 6],...
    'box','off')

print(f, sprintf('%s/JOG-19-0112.Figure2.pdf',OutputFolder), '-dpdf','-bestfit')
    if strcmp(vis,'off')
        close(f)
    end