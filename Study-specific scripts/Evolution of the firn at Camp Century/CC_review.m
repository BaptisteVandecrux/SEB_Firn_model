LRin_org = ncread('C:/Data_save/CC/GITS_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10_save/surf_CC_RACMO.nc','LRin');
melt_mweq = ncread('C:/Data_save/CC/GITS_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10_save/surf_CC_RACMO.nc','melt_mweq');

figure
subplot(1,2,1)
hold on
scatter(LRin_org, data_surf{1}.LRin(1:length(LRin_org)),'.k')
plot([80 350],[80 350],'r','LineWidth',2)
axis tight square; box on; grid on
xlabel('LRin adjusted using THU_U (W/m-2)')
ylabel('LRin adjusted using CEN (W/m-2)')

subplot(1,2,2)
hold on
scatter(melt_mweq*1000, data_surf{1}.melt_mweq(1:length(meltflux))*1000,'.k')
plot([0 1.30],[0 1.30],'r','LineWidth',2)
axis tight square; box on; grid on
xlabel('Melt adjusted using THU_U (mm)')
ylabel('Melt adjusted using CEN (mm)')

%%
% close all
figure
% subplot(2,1,1)
% hold on
% plot(data_AWS_GITS.time,data_AWS_GITS.AirTemperature1C)
% plot(data_AWS_GITS.time,data_AWS_GITS.AirTemperature2C)
% subplot(2,1,2)
hold on
scatter(data_AWS_GITS.AirTemperature2C,data_AWS_GITS.AirTemperature1C,'.k')
plot([-40 5],[-40 5],'r','LineWidth',2)
xlim([-40 5])
ylim([-40 5])
axis square; box on; grid on
xlabel('Air temperature at the lower level (degC)')
ylabel('Air temperature at the upper level (degC)')
text1 = sprintf('ME = %0.1f degC',...
                nanmean( data_AWS_GITS.AirTemperature2C-...
                data_AWS_GITS.AirTemperature1C));
title(text1)

% figure
% subplot(2,1,1)
% hold on
% plot(data_AWS_GITS.time,data_AWS_GITS.AirTemperature3C)
% plot(data_AWS_GITS.time,data_AWS_GITS.AirTemperature4C)
% subplot(2,1,2)
% hold on
% scatter(data_AWS_GITS.AirTemperature3C,data_AWS_GITS.AirTemperature4C,'.k')
% plot([-40 5],[-40 5],'r','LineWidth',2)
% xlim([-40 5])
% ylim([-40 5])
% axis square; box on; grid on
% xlabel('Air temperature at the lower level (degC)')
% ylabel('Air temperature at the upper level (degC)')
% text1 = sprintf('ME = %0.1f degC',...
%                 nanmean( data_AWS_GITS.AirTemperature3C-...
%                 data_AWS_GITS.AirTemperature4C));
% title(text1)
