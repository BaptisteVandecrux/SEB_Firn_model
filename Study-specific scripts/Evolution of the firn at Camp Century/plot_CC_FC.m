clear all
close all
set(0,'defaultfigurepaperunits','centimeters');
   set(0,'DefaultAxesFontSize',16)
   set(0,'defaultfigurecolor','w');
set(0,'defaultfigureinverthardcopy','off');
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[29.7 25]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 25]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 25]-0.5]);
addpath(genpath('.\lib'))
addpath(genpath('Input'),genpath('Output'))

WorkingFolder = './Output';
mkdir(sprintf('%s/Plots',WorkingFolder));
OutputFolder = sprintf('%s/Plots',WorkingFolder);
list = dir(WorkingFolder);
folder_list = {list.name};
folder_list = {folder_list{3:end}};
folder_list(find(strcmp(folder_list,'Plots')))= [];
folder_list(find(strcmp(folder_list,'distrib')))= [];
path_list={};
station={};
 count = 0;
 for i = [1 2 4 5 6 7 3 8 9] % 1:length(folder_list)
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

vis = 'on';
load('CC_FAC_20_save.mat')

%% Plotting CC
f=figure('Visible',vis,'outerposition',[1 -1  20 15]);
%  ha = tight_subplot(3,3,[0.02 0.05],[0.15 0.02],[0.07 0.1]);
 ha = tight_subplot(1,2,0.02,0.2,[0.1 0.18]);
summary_CC = NaN(9,9);
ps = @(m,rho) max(0,m.*(1./rho -1/917));

col = brewermap(9,'paired');
for ii = 1:length(path_list)
   
    set(f,'CurrentAxes',ha(1))
    hold on
    plot(time_mod{ii}, CC_20{ii}/334,'LineWidth',2,'Color',col(ii,:))
    lm = fitlm(time_mod{ii}, CC_20{ii});
    summary_CC(1,ii) = mean(CC_20{ii});
    summary_CC(2,ii) = lm.Coefficients.Estimate(2)*10;
    summary_CC(3,ii) = coefTest(lm);

        set(f,'CurrentAxes',ha(2))
    hold on
    plot(time_mod{ii}, FC20{ii},'LineWidth',2,'Color',col(ii,:))
    disp(station{ii})
    disp((FC20{ii}(end)-FC20{ii}(1))/FC20{ii}(1))
end
    set(f,'CurrentAxes',ha(1))
axis tight
ylabel('Refreezing capacity in top 20 m (mm)')
xlabel('Year')
box on
ha(1).XTick = datenum(1994:2:2016,1,1);
ha(1).XMinorTick = 'on';
datetick('x','yyyy','keepticks','keeplimits')
ha(1).XTickLabelRotation = 45;

    set(f,'CurrentAxes',ha(2))
axis tight
ylabel('Firn capacity in top 20 m (mm)','Interpreter','tex')
box on
xlabel('Year')
ha(2).XTick = datenum(1994:2:2016,1,1);
ha(2).XMinorTick = 'on';
ha(2).XTickLabelRotation = 45;
datetick('x','yyyy','keepticks','keeplimits')
ha(2).YAxisLocation = 'right';
    ha(1).TickLength = 2*ha(1).TickLength;
    ha(2).TickLength = 2*ha(2).TickLength;

legendflex(station,'ref',gcf,'anchor',{'n','n'},'nrow',3)

print(f, sprintf('%s/ColdContent_20',OutputFolder), '-dpng')
