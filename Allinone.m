%% NetCDF file reader

ncfile_firn_density = ('./Output/Kan_U_BIGiceslab_100layers/rho_bin_1.nc');
ncfile_ice_content = ('./Output/Kan_U_BIGiceslab_100layers/snic_bin_1.nc');
ncfile_snow_content = ('./Output/Kan_U_BIGiceslab_100layers/snowc_bin_1.nc');
rho = ncread(ncfile_firn_density,'rho') ;
Depth = ncread(ncfile_firn_density,'Depth') ;
ice = ncread(ncfile_ice_content,'snic');
snow = ncread(ncfile_snow_content,'snowc');
DP_2012 = readtable('RetMIP_density_KAN-U 2012.csv');


%% Calculation for bulk density. 

rho_all = (snow + ice)./(snow./rho + ice/917);

%% Variable generator

depth_layers_final = Depth(:,end);
depth_layers_init = Depth(:,1);
rho_sim_init = rho_all(:,1);
rho_sim_final = rho_all(:,end);


%% 2012 DATA

AVGD_per_metre_2012 = table2array(readtable('2012AVGD_per_metre.txt'));
AVGD_per_halfmetre_2012 = table2array(readtable('2012AVGD_per_halfmetre.txt'));
AVGD_per_2metres_2012 = table2array(readtable('2012AVGD_per_2metres.txt'));
AVGD_per_5metres_2012 = table2array(readtable('2012AVGD_per_5metres.txt'));
AVGD_per_10cm_2012 = table2array(readtable('2012AVGD_per_10cm.txt'));

%% Interpolator final time-step sim values - Take top 20m of firn, round to nearest cm and interpolate

Top20_Depth_final = depth_layers_final(1:82);
Rounded_Depth_final= round(Top20_Depth_final,2);
y = linspace(min(Rounded_Depth_final),max(Rounded_Depth_final),200);
rounded_again_final = round(y,2);
rounded_Depthfinal = rounded_again_final';

Top20_Densities_final = rho_sim_final(1:82);

x_Depths_Final = Rounded_Depth_final;
y_Densities_Final = Top20_Densities_final;
xi_interpolated_depths_final = rounded_Depthfinal;
yi_interpolated_densities_final = interp1(x_Depths_Final,y_Densities_Final,xi_interpolated_depths_final,'next');

%% Interpolator final time-step sim values - Take top 20m of firn, round to nearest cm and interpolate

Top20_Depth_init = depth_layers_init(1:82);
Rounded_Depth_init= round(Top20_Depth_init,2);
y = linspace(min(Rounded_Depth_init),max(Rounded_Depth_init),200);
rounded_again_init = round(y,2);
rounded_Depthinit = rounded_again_init';

Top20_Densities_init = rho_sim_init(1:82);

x_Depths_init = Rounded_Depth_init;
y_Densities_init = Top20_Densities_init;
xi_interpolated_depths_init = rounded_Depthinit;
yi_interpolated_densities_init = interp1(x_Depths_init,y_Densities_init,xi_interpolated_depths_init,'next');


%% Interpolation plot checker final time step sim

subplot(1,2,1)
plot(Top20_Densities_final,Top20_Depth_final)
set(gca,'Ydir','reverse')
title('Simulated density profile')
xlabel('Density kg m^{⁻3}')
ylabel('Depth(m)')
subplot(1,2,2)
plot(y_Densities_Final,x_Depths_Final,'o',yi_interpolated_densities_final,xi_interpolated_depths_final)
set(gca,'Ydir','reverse')
title('Interpolated density profile of simulated data')
xlabel('Density kg m^{⁻3}')
ylabel('Depth(m)')
legend('Interpolated data points')
saveas(gcf,'Simulated density profile, Top 20m, final timestep & interpolated data','png')

%% Interpolation plot checker initial time step sim

% set(findall(fig, '-property', 'fontsize'), 'fontsize', NewFontSize)
% fig1 = gcf;
subplot(1,2,1)
plot(Top20_Densities_init,Top20_Depth_init)
set(gca,'Ydir','reverse')
title('Simulated density profile')
xlabel('Density kg m^{⁻3}')
ylabel('Depth(m)')
subplot(1,2,2)
plot(y_Densities_init,x_Depths_init,'o',yi_interpolated_densities_init,xi_interpolated_depths_init)
set(gca,'Ydir','reverse')
title('Interpolated density profile of simulated data')
xlabel('Density kg m^{⁻3}')
ylabel('Depth(m)')
legend('Interpolated data points')
% saveas(gcf,'Simulated density profile, Top 20m, initial timestep & interpolated data','png')

%% Discretizing and RMS
% simulation contains the variables depth_layers_final, ind (resolution),
% rho_sim_final, rho_sim_init, 'res' sets resolution of avg dens per 'r'

indices = discretize(xi_interpolated_depths_final,res);
good_ind = 1:(find(isnan(indices),1,'first')-1);
rho_simulated_finalDP = accumarray(indices(good_ind),yi_interpolated_densities_final(good_ind),[],@mean);

%% This does all of what happens above but in less code via the function.
% TO DO - 'r' can be a string in the file name i.e simulated_finalrho_AVGDP'r'M
% sometimes there can be a zero entry as the depth values might skip some
% metres in the deeper parts of the firn. 

r = 1;
res = 0:r:20;

simulated_finalrho_AVGDP1M = ConvertDepth(xi_interpolated_depths_final,yi_interpolated_densities_final,res);
simulated_initialrho_AVGDP1M = ConvertDepth(xi_interpolated_depths_init,yi_interpolated_densities_init,res);
A2M = numel(simulated_initialrho_AVGDP1M);
B2M = numel(simulated_finalrho_AVGDP1M);


%% RMS difference between final and initial simulated profiles

RMSD_per10cm_sims = mean((simulated_finalrho_AVGDP10cm-simulated_initialrho_AVGDP10cm).^2);
RMSD_perhalfmeter_sims = mean((simulated_finalrho_AVGDPHM-simulated_initialrho_AVGDPHM).^2);
RMSD_permeter_sims = mean((simulated_finalrho_AVGDP1M-simulated_initialrho_AVGDPM).^2);
RMSD_per2meters_sims = mean((simulated_finalrho_AVGDP2M-simulated_initialrho_AVGDP2M).^2);
RMSD_per5meters_sims = mean((simulated_finalrho_AVGDP5M-simulated_initialrho_AVGDP5M).^2);

%% RMS difference between final simulated profile and 2012 profile at KANU

RMSD_per_meter_2012 = mean((simulated_finalrho_AVGDP1M-AVGD_per_metre_2012).^2);

%%
figure
hold on
stairs(depth_in1, rho_in1)
stairs(depth_in2, rho_in2)


%% Plot and data

% initial and final density profiles
hold on
plot(rho_sim_final,depth_layers_final)
plot(rho_sim_init,depth_layers_init)
set(gca,'Ydir','reverse')
hold on
title('Initial and final simulated density profiles')
legend('final profile', 'initial profile','Location','west')
xlabel('density (kg/m^{-3})')
ylabel('depth(m)')
%saveas(gcf,'Initial and final simulated density profiles','png')


%% Plot distribution of layers per depth

hold on
histogram(depth_layers_final,101);
histogram(depth_layers_init,101);
T20layers_fin = depth_layers_final<=20;
T20layers_init = depth_layers_init <=20;
T20Layerscount_fin = sum(T20layers_fin);
T20Layerscount_init = sum(T20layers_init);
title('Distribution of layers in top 20m of firn')
legend('simulated final', 'simulated initial')
str = {'Number of layers in final time step =',T20Layerscount_fin, 'Number of layers in initial time step=',T20Layerscount_init};
text(20,5,str,'FontSize',10)
%saveas(gcf,'Layer distribution histograms','png')



% %% Plot 2012 profile  against simulated profiles 
% %% NEEDS TO BE INTERPOLOATED FOR SIMULATED PLOTS
% 
% hold on
% plot(DP_2012{1:200,2},DP_2012{1:200,1})
% %plot(rho_sim_final,depth_layers_final)
% %plot(rho_sim_init,depth_layers_init)
% %ylim([0 62])
% xlabel('Density(kg m^{⁻3})');
% ylabel('Depth(m)');
% title('KANU 2012 profile and simulated initial and final profiles')
% set(gca, 'YDir','reverse') 
% legend('KANU 2012', 'initial sim', 'final sim')

