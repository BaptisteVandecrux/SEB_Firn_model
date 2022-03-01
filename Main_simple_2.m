% Main script for running the surface - subsurface model
% Here you can define which year you want to compute, define your
% param{kk}eters, plot output variables, evaluate model performance....
%
% Author: Baptiste Vandecrux (bava@byg.dtu.dk)
% ========================================================================
clearvars
close all
clc
addpath(genpath('.\lib'))
addpath(genpath('Input'),genpath('Output'))

% Running a single file
NumLayer = 50;
param.z_max = 50;
param.dz_ice = param.z_max/NumLayer;
param.verbose = 1;
param.lim_new_lay = param.z_max/NumLayer/50;
param.shallow_firn_lim = 5;

param.retmip = 1;
param.vis = 'off';

model_version = 'imau_antarctica';
station_list = {'IMAU_aws4'};
for i =1:length(station_list)
    param.station = station_list{i};
    param.InputAWSFile = ['./Input/', station_list{i}, '_high-res_meteo.csv'];
%     param.InputAWSFile = 'Input\Weather data\data_KAN_U_2012-2016.txt';
    [RunName, c] = HHsubsurf(param);
    disp(RunName)
    OutputFolder = ['Output/', RunName];
    vis = 'on';
    ylimit = 10;
    load(strcat(OutputFolder,'/run_param.mat'))
    c.OutputFolder = OutputFolder;
    namefile = sprintf('%s/%s_surface.nc',OutputFolder,c.station);
    finfo = ncinfo(namefile);
    names={finfo.Variables.Name};
    for i= 1:size(finfo.Variables,2)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end

    % extract subsurface variables
    varname= {'compaction' 'dgrain' 'rfrz' 'rho' 'slwc' 'snowc' 'snic' 'T_ice'};
    for i =1:length(varname)
        namefile = sprintf('%s/%s_%s.nc',...
            OutputFolder,c.station,varname{i});

        eval(sprintf('%s = ncread(''%s'',''%s'');', varname{i}, namefile,varname{i}));
    end
    fprintf('\nData extracted from nc files.\n');

    thickness_weq = snowc + snic +slwc;

    % Update BV2017: (not used anymore) Liquid water does not participate to the actual thickness
    % of a layer. Ice does not participate as long as there is enough pore
    % space in the snow to accomodate it.
    pore_space = snowc .* c.rho_water.*( 1./rho - 1/c.rho_ice);
    excess_ice = max(0, snic * c.rho_water / c.rho_ice - pore_space);
    thickness_act = snowc.*(c.rho_water./rho) + excess_ice;

    depth_act=zeros(size(thickness_act));
    depth_weq=zeros(size(thickness_weq));
    for j=1:length(thickness_act(1,:))
            depth_act(1:end,j)=cumsum(thickness_act(1:end,j));
            depth_weq(1:end,j)=cumsum(thickness_weq(1:end,j));
    end

    rho_all= (snowc + snic)./...
                (snowc./rho + snic./c.rho_ice);
    lwc = slwc(:,:) ./ thickness_act(:,:);

    time_mod = datenum(Year,1,Day,Hour,0,0);
    TT = ones(c.jpgrnd,1) * time_mod';

    H_surf_old = H_surf;
    H_surf = depth_act(end,:)'-depth_act(end,1); %+snowbkt*1000/315;

    for i = 1:length(H_surf)-1
        if (H_surf(i+1)-H_surf(i))> c.new_bottom_lay-1
            H_surf(i+1:end) = H_surf(i+1:end) - c.new_bottom_lay*c.rho_water/c.rho_ice;
        end
    end    

    if ~isfield(c,'verbose')
        c.verbose = 1;
    end
    % extracting observed

    [data_AWS,Tsurf_obs, pres, T1, T2, z_T1, z_T2, ...
        o_T1, o_T2,RH1, RH2, z_RH1, z_RH2, ...
        o_RH1, o_RH2,WS1, WS2, ~, z_WS2, ...
        o_WS1, o_WS2, SRin, SRout, LRin, LRout, time_yr, year, day, hour, c] = PrepareForRetMIP(c);
        T_ice_obs = [];
        depth_thermistor = [];


    time_obs = datenum(year,1,day,hour,0,0);
    depth_obs = depth_thermistor';
    depth_obs(depth_obs==0) = NaN;
    depth_act_save = depth_act;
    depth_act = vertcat(zeros(size(depth_act(1,:))), depth_act);
    depth_weq = vertcat(zeros(size(depth_weq(1,:))), depth_weq);

    for i = 1:size(depth_act,1)
        depth_act(i,:) = depth_act(i,:) - H_surf' + H_surf(1);
    end

    disp('Plotting water content')
    f = figure('Visible', vis);
    ha = tight_subplot(2,1,0.07,0.07,0.07);
    set(ha(2),'Visible','off')
    set(f,'CurrentAxes',ha(1))
    col = PlotTemp(TT,depth_act,slwc*1000,...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','cool',...
        'XLabel','Time',...
        'YLabel','Depth (m)',...
        'CLabel','\theta_{liq} (mm w.eq.)',...
        'Range', 0:0.1:1);
    if  strcmp(c.station, 'GITS')
        ylim([-5 ylimit])
    else
        plot(time_mod,-H_surf+H_surf(1),'LineWidth',2)
        ylim([min([-H_surf+H_surf(1)])-1 ylimit])
    end

    col.TickLength=col.TickLength*5;
    temp = col.YTickLabel;
    for i = 1:length(temp)
        if i/2==floor(i/2)
            temp(i,:)=' ';
        end
    end

    set(gca,'Color',[0.95 0.95 0.95]);
    print(f, sprintf('%s/slwc',OutputFolder), '-djpeg')
    if strcmp(vis,'off')
        close(f);
    end
end