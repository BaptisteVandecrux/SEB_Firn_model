% Model run analysis.
% Here you can plot output variables, evaluate model performance, etc...
%
% Author: Baptiste Vandecrux (bava@byg.dtu.dk)
% ========================================================================
clearvars
close all
clc
addpath(genpath('lib'))
addpath(genpath('Input'),'Output')

set(0,'defaultfigurepaperunits','centimeters');
set(0,'DefaultAxesFontSize',15)
set(0,'defaultfigurecolor','w');
set(0,'defaultfigureinverthardcopy','off');
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperpositionmode','auto');
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
set(0,'defaultfigurepapersize',[29.7 16]);

PlotThings = 'yes'; % enable/disable all plots
vis = 'off';        % make the plot visible or not (only plotted in files)

% ========= Upload model results and calculating depth scales ============
% Read surface variables
    OutputFolder = uigetdir('./Output');
%     OutputFolder = '.\Output\GITS_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10';

%     folderlist = { ...
%         ...'CP1_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10',...
%         'DYE-2_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10',...
%        	'Summit_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10',...
%         'NASA-SE_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10'};
% for i=1:length(folderlist)
%     OutputFolder = sprintf('./Output/%s',folderlist{i})
%     close all
    

    
    % extract run parameters
    load(strcat(OutputFolder,'/run_param.mat'))
    c.OutputFolder = OutputFolder;
    % extract surface variables
    switch c.station
        case 'DYE-2'
            namefile = sprintf('%s/%s_surf-bin-1.nc',OutputFolder,c.station);
        otherwise
            namefile = sprintf('%s/surf-bin-%i.nc',OutputFolder,1);
    end
    finfo = ncinfo(namefile);
    names={finfo.Variables.Name};
    for i= 1:size(finfo.Variables,2)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end

    % extract subsurface variables
    varname= {'compaction' 'dgrain' 'rfrz' 'rho' 'slwc' 'snowc' 'snic' 'T_ice'};
    for i =1:length(varname)
        switch c.station
            case 'DYE-2'
                namefile = sprintf('%s/%s_%s_bin_1.nc',...
                    OutputFolder,c.station,varname{i});
        end
        try eval(sprintf('%s = ncread(''%s'',''%s'');', varname{i}, namefile,varname{i}));
        catch me
            namefile = sprintf('%s/%s_bin_%i.nc',OutputFolder,varname{i},1);
            eval(sprintf('%s = ncread(''%s'',''%s'');', varname{i}, namefile,varname{i}));
        end
    end
    fprintf('\nData extracted from nc files.\n');


thickness_weq = snowc + snic +slwc;

% Update BV2017: (not used anymore) Liquid water does not participate to the actual thickness
% of a layer. Ice does not participate as long as there is enough pore
% space in the snow to accomodate it.
pore_space = snowc .* c.rho_water.*( 1./rho - 1/c.rho_ice);
excess_ice = max(0, snic * c.rho_water / c.rho_ice - pore_space);
thickness_act = snowc.*(c.rho_water./rho) + excess_ice;

% thickness_act = snowc.*(c.rho_water./rho) + snic.*(c.rho_water./c.rho_ice);

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

[time_yr, year, day, hour, pres,...
    T1, T, z_T1, H_instr_temp, o_T1,o_T2, ...
    RH1, RH, z_RH1, H_instr_hum, o_RH1, o_RH2, ...
    WS1, WS, z_WS1, H_instr_wind, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_AWS, c] = ...
    ExtractAWSData(c);

time_obs = datenum(year,1,day,hour,0,0);
depth_obs = depth_thermistor';
depth_obs(depth_obs==0) = NaN;

T = T' + c.T_0;
T(isnan(depth_obs)) = NaN;
TT_obs= repmat(time_obs',size(depth_obs,1),1);
% depth scale
depth_act_save = depth_act;
depth_act = vertcat(zeros(size(depth_act(1,:))), depth_act);
depth_weq = vertcat(zeros(size(depth_weq(1,:))), depth_weq);

for i = 1:size(depth_act,1)
    depth_act(i,:) = depth_act(i,:) - H_surf' + H_surf(1);
end

% figure
% plot(depth_act_save(:,1),rho_all(:,1))
% hold on
% plot(depth_act_save(:,end),rho_all(:,end))
% title(c.station)
% 
% depth = [0:0.1:depth_act_save(end,end)]';
% density = interp1([0; depth_act_save(:,end)],[315; rho_all(:,end)],depth);
% ice_perc = NaN*density;
% 
% filename = sprintf('./Input/Initial state/DensityProfile_%s_%i.csv',c.station,1998);
%     M  = [depth, density];
%     M_table = array2table(M,'VariableName', {'depth_m', 'density_kgm3'});
% 
%     writetable(M_table,filename,'Delimiter',';')

%% ===================== Outputing text files ======================
% 2 m temperature humidity and 10 m wind speed
disp('Writing 2m T, 2m RH and 10m WS')
time = time_mod;

M = table(time,time_yr, ...
    theta_2m-c.T_0,data_AWS.AirTemperature2C_Origin, ...
    RH_2m_wrtw,...
    data_AWS.RelativeHumidity2_Origin,...
    ws_10m, data_AWS.WindSpeed2ms_Origin);
M.Properties.VariableNames = {'time','time_yr','T_2m_degC','T_Origin', 'RH_2m_perc', ...
    'RH_Origin','WS_10m_ms','WS_Origin'};
writetable(M,sprintf('%s/%s_T_RH_WS_standard.csv',OutputFolder,c.station),'Delimiter',';');
PlotWeather(M, ...
    'VarList',{'T_2m_degC','RH_2m_perc', 'WS_10m_ms'},...
    'LabelList', {'Air temperature \newline          (^oC)',...
    '    Relative \newline humidity (%)',...
    'Wind speed \newline      (m/s)'},...
    'Origins', {'T_Origin','RH_Origin','WS_Origin'},...
    'OutputFolder',OutputFolder,'vis',vis);

% SEB variables 
disp('Writing hourly SEB')
SEB_hour = table(time,time_yr, SHF, LHF, SRin - SRout,LRin - LRout_mdl,...
    rainHF,GF, meltflux,meltflux.*c.dt_obs./c.dev./c.L_fus./c.rho_water, Tsurf,...
    'VariableNames',{'time', 'time_yr','SHF_Wm2','LHF_Wm2','SRnet_Wm2',...
    'LRnet_Wm2','RainHF_Wm2','GF_Wm2','MeltEnergy_Wm2',...
    'Melt_mweq','Tsurf_degC'});
writetable(SEB_hour, sprintf('%s/SEB_hour.txt', OutputFolder))

disp('Writing JJA SEB')
SEB_hour_MJ = SEB_hour;
for i = 3:9
    SEB_hour_MJ.(SEB_hour_MJ.Properties.VariableNames{i}) = ...
            SEB_hour_MJ.(SEB_hour_MJ.Properties.VariableNames{i}).*3600/1000000;
end
SEB_JJA = AvgTableJJA(SEB_hour_MJ,@sum);
writetable(SEB_JJA, sprintf('%s/SEB_JJA.txt', OutputFolder))

disp('Writing hourly, monthly, yearly, seasonal SMB')
% runoffhour = -[runoff(1); runoff(2:end)-runoff(1:end-1)] ;
SMB = snowfall+rainfall+sublimation-runoff;

SMB_hour = table(time, time_yr, snowfall, rainfall, sublimation, runoff,SMB);
SMB_month = AvgTable(SMB_hour,'monthly');
writetable(SMB_hour, sprintf('%s/SMB_hour.txt', OutputFolder));
writetable(SMB_month, sprintf('%s/SMB_month.txt', OutputFolder));

if length(time_obs)>365*24
    SMB_year = AvgTable(SMB_hour,'yearly');
    SMB_season = AvgTable(SMB_hour,'seasonaly');
    SMB_wateryear = AvgTable(SMB_hour,'water-yearly');
    writetable(SMB_year, sprintf('%s/SMB_year.txt', OutputFolder));
    writetable(SMB_season, sprintf('%s/SMB_season.txt', OutputFolder));
    writetable(SMB_wateryear, sprintf('%s/SMB_wateryear.txt', OutputFolder))
end

 %% ==================== Specific Plots =====================================
if strcmp(PlotThings, 'yes')
    if strcmp(c.station,'NUK_K')
        PlotStuff_NUK_K(time_mod, Tsurf_obs,Tsurf, T_ice, ...
            depth_act_save, rho_all, Surface_Height,SMB_wateryear, OutputFolder,vis,c)
    end
end

 %% ==================== Standard Plots =====================================
if strcmp(PlotThings, 'yes')
    T_subsurf_mod = PlotStuff(Tsurf_obs, Tsurf, TT, depth_act,depth_act_save, depth_weq, ...
        T_ice, rho_all, rho, snowc, snic, slwc, rfrz, time_mod, Surface_Height, ...
        snowfall, rainfall, runoff,  H_surf,  H_comp, SHF, LHF, ...
        SRin, SRout, LRin, LRout_mdl, rainHF, meltflux, TT_obs, depth_obs, ...
        T, sublimation, OutputFolder, compaction,...
        thickness_act, thickness_weq, T, GF,dgrain, data_AWS, vis, c);
end

%% ================= Model Performance Analysis ===========================
fileID = fopen(sprintf('%s/log.txt',OutputFolder),'a+');
names = fieldnames(c);
fprintf(fileID, 'List of parameters:\n');

for i=1:length(names)
    if length(c.(names{i}))<100
        if ischar(c.(names{i}))
            fprintf(fileID, '%s\t=\t',names{i});
            fprintf(fileID, '%s   ',c.(names{i}));
            fprintf(fileID, '\n');
        else
            fprintf(fileID, '%s\t=\t',names{i});
            fprintf(fileID, '%e   ',c.(names{i}));
            fprintf(fileID, '\n');
        end
    end
end

fprintf(fileID, '\n');
[RMSE_tot] = ModelPerformance('RMSE', time_mod,time_obs, Tsurf, ...
    LRout, H_surf, depth_obs, Surface_Height, T_subsurf_mod, T, ...
    fileID, c);
[MSE_tot] = ModelPerformance('MSE', time_mod,time_obs, Tsurf, ...
    LRout, H_surf, depth_obs, Surface_Height, T_subsurf_mod, T, ...
    fileID, c);
[SS_tot] = ModelPerformance('SS', time_mod,time_obs, Tsurf, ...
    LRout, H_surf, depth_obs, Surface_Height, T_subsurf_mod, T, ...
    fileID, c);
[NashSut_tot] = ModelPerformance('NashSut', time_mod,time_obs, Tsurf, ...
    LRout, H_surf, depth_obs, Surface_Height, T_subsurf_mod, T, ...
    fileID, c);

disp('End analysis.')
fclose(fileID);
diary off
% end

