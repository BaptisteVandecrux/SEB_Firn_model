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

model_list = {'RACMO',  'CanESM2_rcp26', 'CanESM2_rcp45', 'CanESM2_rcp85'};
model_list2 = {'RACMO';  'CanESM RCP2.6'; 'CanESM RCP4.5'; 'CanESM RCP8.5'};
path_list = {'C:/Data_save/CC/GITS_RACMO_GEUS_model_output',...
 'C:/Data_save/CC/GITS_CanESM2_RCP26_GEUS_model_output',...
 'C:/Data_save/CC/GITS_CanESM2_RCP45_GEUS_model_output',...
 'C:/Data_save/CC/GITS_CanESM2_RCP85_GEUS_model_output'};

OutputFolder = './Output/Camp Century';
mkdir(OutputFolder)
load ../matlab_functions/Core_all
abc = char(97:122) ;
ABC = char(65:90) ;

%%  Upload model results and calculating depth scales 
step = 1;

for ii = 1:4
    disp(model_list{ii})
    tic
    % extract run parameters
    load(strcat(path_list{ii},'/run_param.mat'));
    c_all{ii} = c;
    disp(c.InputAWSFile)
    % extract surface variables
    namefile = sprintf('%s/CC_%s_surface.nc',path_list{ii},model_list{ii});
    names =  {'time' ...
        'LRin' 'LRout_mdl' 'SHF' 'LHF' 'GF' 'rainHF' 'meltflux'...
        'H_surf' 'SMB_mweq' 'melt_mweq' 'sublimation_mweq' ...
        'H_comp' 'runoff' 'snowthick'...
        'snowfall' 'rainfall' 'SRin','SRout', 'Tsurf'...
        'theta_2m','RH_2m_wrtw','ws_10m'};
    
    data_surf{ii}= table();
    for i= 1:length(names)
        eval(sprintf('data_surf{ii}.%s = ncread(''%s'',''%s'',1,Inf,%i);',...
            char(names{i}), namefile,char(names{i}), step));
    end

    % extract subsurface variables
    varname= {'compaction'  'rho_firn_only' 'slwc' 'snowc' 'snic' 'T_ice'};
    data_subsurf{ii} = table();

    for i =1:length(varname)
        namefile = sprintf('%s/CC_%s_%s.nc',path_list{ii},model_list{ii},...
            varname{i});
        if i==2
            varname{i}='rho';
        end
        eval(sprintf('data_subsurf{ii}.%s = ncread(''%s'',''%s'',[1 1], Inf*[1 1],[1 %i]);',...
            varname{i}, namefile,varname{i},step));
    end
    fprintf('\nData extracted from nc files.\n');
    data_surf{ii}.snowfall = data_surf{ii}.snowfall*step;
    data_surf{ii}.rainfall = data_surf{ii}.rainfall*step;
    data_surf{ii}.sublimation_mweq = data_surf{ii}.sublimation_mweq*step;
    data_surf{ii}.runoff = data_surf{ii}.runoff*step;
    data_subsurf{ii}.compaction = data_subsurf{ii}.compaction*step;
    
    data_subsurf{ii}.thickness_weq = ...
        data_subsurf{ii}.snowc + data_subsurf{ii}.snic + data_subsurf{ii}.slwc;

    % Update BV2017: (not used anymore) Liquid water does not participate to the actual thickness
    % of a layer. Ice does not participate as long as there is enough pore
    % space in the snow to accomodate it.
    pore_space = data_subsurf{ii}.snowc .* c.rho_water.*( 1./data_subsurf{ii}.rho - 1/c.rho_ice);
    excess_ice = max(0, data_subsurf{ii}.snic * c.rho_water / c.rho_ice - pore_space);
    data_subsurf{ii}.thickness_act = data_subsurf{ii}.snowc.*(c.rho_water./data_subsurf{ii}.rho) + excess_ice;

    % data_subsurf{ii}.thickness_act = snowc.*(c.rho_water./rho) + snic.*(c.rho_water./c.rho_ice);

    length_time{ii} = length(data_surf{ii}.time);
    num_layer{ii} = size(data_subsurf{ii}.thickness_weq,1);
    
    data_subsurf{ii}.depth_act = zeros(num_layer{ii}, length_time{ii});
    data_subsurf{ii}.depth_act=...
        cumsum(data_subsurf{ii}.thickness_act,1);

    data_subsurf{ii}.rho_all= (data_subsurf{ii}.snowc + data_subsurf{ii}.snic)./...
                (data_subsurf{ii}.snowc./data_subsurf{ii}.rho + data_subsurf{ii}.snic./c.rho_ice);
    data_subsurf{ii}.lwc = data_subsurf{ii}.slwc(:,:) ./ data_subsurf{ii}.thickness_act(:,:);

    data_surf{ii}.time = ...
        data_surf{ii}.time + datenum(1900,1,1)-1;
    data_surf{ii}.time =  datenum(0,1,data_surf{ii}.time(1),...
        0:(length(data_surf{ii}.time)-1),0,0)';
    data_subsurf{ii}.TT = ones(c.jpgrnd,1) * data_surf{ii}.time';

%     data_surf{ii}.H_surf_old = data_surf{ii}.data_surf{ii}.H_surf;
    data_surf{ii}.H_surf = data_subsurf{ii}.depth_act(end,:)'-data_subsurf{ii}.depth_act(end,1); %+snowbkt*1000/315;

    for i = 1:length_time{ii}-1
        if (data_surf{ii}.H_surf(i+1)-data_surf{ii}.H_surf(i))> c.new_bottom_lay-1
            data_surf{ii}.H_surf(i+1:end) = ...
                data_surf{ii}.H_surf(i+1:end) - ...
                c.new_bottom_lay*c.rho_water/c.rho_ice;
        end
    end  
    % depth scale
    data_subsurf{ii}.depth_act_m = data_subsurf{ii}.depth_act- ...
        repmat(data_surf{ii}.H_surf',num_layer{ii},1);
    
    data_surf{ii}.theta_2m = data_surf{ii}.theta_2m-c_all{1}.T_0;
    data_surf{ii}.SRnet = data_surf{ii}.SRin - data_surf{ii}.SRout;
    data_surf{ii}.LRnet = data_surf{ii}.LRin - data_surf{ii}.LRout_mdl;
    data_surf{ii}.albedo = data_surf{ii}.SRout ./ data_surf{ii}.SRin;
    DV = datevec(data_surf{ii}.time);
    data_surf{ii}.albedo(~ismember(DV(:,2),[6 7 8])) = NaN;
    data_surf{ii}.albedo(data_surf{ii}.albedo>1) = 1;
    data_surf{ii}.albedo(data_surf{ii}.albedo<0.4) = 0.4;

    data_surf_yr{ii} = AvgTable(data_surf{ii},'yearly',@nanmean);
    data_surf_yr2{ii} = AvgTable(data_surf{ii},'yearly',@sum);
    data_surf_yr{ii}.snowfall=data_surf_yr2{ii}.snowfall;
    data_surf_yr{ii}.rainfall =data_surf_yr2{ii}.rainfall;
    data_surf_yr{ii}.sublimation_mweq=data_surf_yr2{ii}.sublimation_mweq;
    data_surf_yr{ii}.melt_mweq=data_surf_yr2{ii}.melt_mweq;
    data_surf_yr{ii}.runoff=data_surf_yr2{ii}.runoff;
    data_surf_yr{ii}.SMB_mweq=data_surf_yr2{ii}.SMB_mweq;
    toc
end
clearvars var_name data_surf_yr2 excess_ice c finfo pore_space Tsurf_obs namefile names Surface_Height

%% 1-Validation accumulation
% CC_10 accumulation
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [6, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["Depth", "Density"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
CC10density = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Data\Camp Century\KU\CC_2010_density.txt", opts);
clear opts

opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["depth_m", "year", "Var3"];
opts.SelectedVariableNames = ["depth_m", "year"];
opts.VariableTypes = ["double", "double", "char"];
opts = setvaropts(opts, 3, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 3, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
CC10annual = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Data\Camp Century\KU\CC10_annual.txt", opts);
clear opts

CC10_interp = table();
CC10_interp.depth= [2.95:0.01:max(CC10density.Depth)]';
CC10_interp.density = interp1([1.1; CC10density.Depth],...
    [NaN; CC10density.Density],CC10_interp.depth,'next');

% figure
% stairs0([NaN; CC10density.Density],...
%     [1.1; CC10density.Depth])
% hold on
% plot(CC10_interp.density,CC10_interp.depth)
% for i = 1:length(CC10annual.depth_m)
%     plot([0 800],[1 1]*CC10annual.depth_m(i))
% end
% set(gca,'YDir','reverse')

subs = discretize(CC10_interp.depth,CC10annual.depth_m);
cc10accum= table();
cc10accum.year = [2006:-1:min(CC10annual.year)]';
cc10accum.b_mweq = accumarray(subs,CC10_interp.density*0.01)/1000;
cc10accum = flipud(cc10accum);
% figure
% stairs(cc10accum.year, cc10accum.b_mweq)

writetable(cc10accum,'Input/Extra/Camp Century/CC10_accum.csv','Delimiter',';')

clearvars h
file_list =     {'CC10_accum.csv'...
    'CC77-1_accum.csv',...
    'CC77-2_accum.csv',...
    'GITS-1_accum.csv',...
    'GITS-2_accum.csv'};

core_name =     {'CC10'...
    'CC77-1',...
    'CC77-2',...
    'GITS-1',...
    'GITS-2'};

accum = cell(1,length(file_list));
for i = 1:length(file_list)
    opts = delimitedTextImportOptions("NumVariables", 2);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ";";
    opts.VariableNames = ["year", "b_mweq"];
    opts.VariableTypes = ["double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    accum{i} = readtable(['Input/Extra/Camp Century/' file_list{i}], opts);
    clear opts
    if i>1
        accum{i}.b_mweq=accum{i}.b_mweq;
    end
end

disp([accum{2}.year([182 212]) accum{3}.year([105 135])])
disp([mean(accum{2}.b_mweq(182:212)) mean(accum{3}.b_mweq(105:135))]*1000/917)
% 1943-73 Mean rate(m of ice / year): 0.380

accum_all = table();
accum_all.year = [1960:2006]';
accum_all.b_mweq = NaN(length(accum_all.year),length(file_list));
for i = 1:length(file_list)
     ind_in_core = ismember(accum{i}.year,accum_all.year);
     ind_in_all = ismember(accum_all.year,accum{i}.year);
     accum_all.b_mweq(ind_in_all,i) = accum{i}.b_mweq(ind_in_core);
end

accum_all.min = nanmin(accum_all.b_mweq,[],2);
accum_all.max = nanmax(accum_all.b_mweq,[],2);
accum_all.median = nanmedian(accum_all.b_mweq,2);
% accum_all.mean = accum_all.b_mweq(:,1);
accum_all.mean = nanmean(accum_all.b_mweq,2);

accum_all_dec = table();
accum_all_dec.year = [1960:10:2000]';
subs = discretize(accum_all.year,[accum_all_dec.year; accum_all_dec.year(end)+10]-0.1);
accum_all_dec.b_mweq = accumarray(subs,accum_all.mean,size(accum_all_dec.year),@nanmean);

accum_all_dec = [accum_all_dec; accum_all_dec(end,:)];
accum_all_dec.year(end) = 2010;

accum_all = [accum_all; accum_all(end,:)];
accum_all.year(end) = 2007;

opts = spreadsheetImportOptions("NumVariables", 11);
opts.Sheet = "Accumulation rates";
opts.DataRange = "D4:N4";
opts.VariableNames = ["VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14"];
opts.SelectedVariableNames = ["VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
BuchardtetalS1 = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Data\Accumulation\Buchardt et al. 2012. Greenland decadal averages of accumulation and temperature 1900-2009.xlsx", opts, "UseExcel", false);
Buchardt = table();
Buchardt.year = [1900:10:2010]';
Buchardt.b_mweq = [table2array(BuchardtetalS1)'; BuchardtetalS1.VarName14(end)]*914/1000;
clear opts

f = figure('Visible',vis,'Outerposition',[1 1 18 23]);
time_accum =1960:2010;
hold on 
% col = lines(length(accum)+1);    
for i =1:length(accum)
    stairs([accum{i}.year; accum{i}.year(end)+1], [accum{i}.b_mweq; accum{i}.b_mweq(end)],... 
        'LineWidth',2);
end  
stairs(Buchardt.year,Buchardt.b_mweq,'b','LineWidth',3)
stairs(accum_all_dec.year,accum_all_dec.b_mweq,'--r','LineWidth',3)
axis tight
    xlim(time_accum([1 end]))
 
tmp = num2str(time_accum(1:5:end)');
for i=2:2:size(tmp,1) 
    tmp(i,:)='    ';
end

set (gca,'layer','top',...
   'YMinorTick','on','XMinorTick','on',...
   'XTick',time_accum(1:5:end),...
   'TickLength', [0.01 0.025].*2,...
   'XTickLabel',tmp,...
   'XTickLabelRotation',0,...
   'layer','top')
legend([core_name 'Buchardt et al. (2012)' 'decadal mean from all 5 cores'],'location','northoutside')
ylabel( sprintf('SMB (m w.e.)'),'Interpreter','tex');
xlabel('Year')
box on; grid on
print(f, sprintf('%s/accum_core',OutputFolder), '-djpeg','-r300')

% f = figure('Visible',vis,'Outerposition',[1 1 40 18]);
% 
% subplot(1,3,1)
% lm1 = fitlm(accum_all.b_mweq(:,1),accum_all.b_mweq(:,2));
% scatter(accum_all.b_mweq(:,1),accum_all.b_mweq(:,2),'fill')
% hold on
% scatter(accum_all.b_mweq(:,1),accum_all.b_mweq(:,3),'fill')
% plot([0.2 0.8],lm1.Coefficients.Estimate(1)+lm1.Coefficients.Estimate(2)*[0.2 0.8],'--b')
% plot([0.2 0.8],[0.2 0.8],'k')
% axis tight square; box on
% xlabel('CC_10 accumulation (m w.e.)')
% ylabel('CC_77 accumulation (m w.e.)')
% legend('using CC_77-1','using CC_77-2','interpreter','none','location','northoutside')
% 
% subplot(1,3,2)
% lm1 = fitlm(accum_all.b_mweq(:,2),accum_all.b_mweq(:,4));
% scatter(accum_all.b_mweq(:,2),accum_all.b_mweq(:,4),'fill')
% hold on
% scatter(accum_all.b_mweq(:,3),accum_all.b_mweq(:,5),'fill')
% plot([0.2 0.8],lm1.Coefficients.Estimate(1)+lm1.Coefficients.Estimate(2)*[0.2 0.8],'--b')
% plot([0.2 0.8],[0.2 0.8],'k')
% axis tight square; box on
% xlabel('CC_77 accumulation (m w.e.)')
% ylabel('GITS accumulation (m w.e.)')
% legend('using core 1','using core 2','location','northoutside')
% 
% subplot(1,3,3)
% lm1 = fitlm(accum_all.b_mweq(:,1),accum_all.b_mweq(:,4));
% scatter(accum_all.b_mweq(:,1),accum_all.b_mweq(:,4),'fill')
% hold on
% scatter(accum_all.b_mweq(:,1),accum_all.b_mweq(:,5),'fill')
% plot([0.2 0.8],lm1.Coefficients.Estimate(1)+lm1.Coefficients.Estimate(2)*[0.2 0.8],'--b')
% plot([0.2 0.8],[0.2 0.8],'k')
% axis tight square; box on
% legend('using GITS-1','using GITS-2','location','northoutside')
% xlabel('CC_10 accumulation (m w.e.)')
% ylabel('GITS accumulation (m w.e.)')
% print(f, sprintf('%s/accum_core_2',OutputFolder), '-djpeg','-r300')



f = figure('Visible',vis,'Outerposition',[1 1 18 15]);
ha = tight_subplot(1,1,0.07,[0.2 0.2],0.1);
set(gcf,'CurrentAxes',ha(1))
time_accum =1966:2010;

hold on 
h(5) = plot(NaN,NaN,'w');
sh1 = stairs(accum_all.year,accum_all.max,'w');
sh2 = stairs(accum_all.year,accum_all.min,'w');
x = [sh1.XData(1),repelem(sh1.XData(2:end),2)];
y1 = [repelem(sh1.YData(1:end-1),2),sh1.YData(end)];
y2 = [repelem(sh2.YData(1:end-1),2),sh2.YData(end)];
x2 = [x, fliplr(x)];
inBetween = [y1, fliplr(y2)];
h(6) =  fill(x2, inBetween, RGB('light gray'),'LineStyle','none');
h(7) = stairs(accum_all.year,accum_all.median,'k','LineWidth',3);

disp(mean(accum_all.median(8:48)))

col_mod = lines(4);
for ii = 1:4
        data_yr = AvgTable(data_surf{ii},'yearly',@sum);
 set(gcf,'CurrentAxes',ha(1))
    DV = datevec(data_yr.time);
    years = DV(:,1);
    month = DV(:,2);
    
    data_yr.snowfall = data_yr.snowfall;
    SMB_station = data_yr.snowfall + data_yr.sublimation_mweq;
%     disp(mean(SMB_station(2:42)))

    h(ii) = stairs0(years+month/12,SMB_station);
    h(ii).Color =    col_mod(ii,:);
    h(ii).LineWidth =    2;
    
    accum_all(accum_all.year<years(1),:) = [];
    fprintf('%i) %s, ME (mm w.e.): %0.2f\n',ii,model_list2{ii},...
        1000*nanmean(SMB_station(1:length(accum_all.median)) - accum_all.median));
    fprintf('relative bias (perc): %0.2f\n',...
        100*nanmean(SMB_station(1:length(accum_all.median))-accum_all.median)/...
        nanmean(SMB_station(1:length(accum_all.median))));
    fprintf('mean sublimation (mm w.e.): %0.2f\n',...
        1000*nanmean(data_yr.sublimation_mweq(1:length(accum_all.median))));
end

axis tight
xlim(time_accum([1 end]))
 
tmp = num2str(time_accum(1:5:end)');
for i=2:2:size(tmp,1) 
    tmp(i,:)='    ';
end

set (gca,'layer','top',...
   'YMinorTick','on','XMinorTick','on',...
   'XTick',time_accum(1:5:end),...
   'TickLength', [0.01 0.025].*2,...
   'XTickLabel',tmp,...
   'XTickLabelRotation',0)

ylabel( sprintf('SMB (m w.e.)'),'Interpreter','tex');
xlabel('Year')
box on; grid on
h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'A');
h_text.FontWeight='bold';
h_text.FontSize=16;

legendflex([h ],[model_list2, 'Observations from 5 firn cores:' 'min-max range','median'],'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 0], ...
                       'ncol',2, ...
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
print(f, sprintf('%s/1-Accumulation_val',OutputFolder), '-djpeg','-r300')

%% 2-Surface forcing, 3-SEB, 4-SMB  

 VarList = {'theta_2m','RH_2m_wrtw','ws_10m','SRin','albedo','LRin'};
LabelList =  {{'Air temperature (^oC)'},...
    {'Relative humidity (%)'},...
    {'Wind speed (m/s)'},...
    {'Downward SW', 'radiation (W/m^2)'},...
    {'Surface albedo (-)'},...
    {'Downward LW', 'radiation (W/m^2)'}};
NameFile = sprintf('%s/_3-SurfaceForcing',OutputFolder);
PlotWeather_CC(data_surf_yr,vis,VarList,LabelList,NameFile,model_list2);
   
VarList = {'SHF','LHF','GF','rainHF','SRnet','LRnet','Tsurf','melt_mweq'};
LabelList =  { {'Mean sensible', 'heat flux (W / m^2)'},...
    {'Mean latent', 'heat flux (W / m^2)'},...
{'Mean conductive', 'heat flux (W / m^2)'},...
{'Mean rain', 'heat flux (W / m^2)'},...
    {'Mean net shortwave', 'radiation (W / m^2)'},...
    {'Mean net longwave', 'radiation (W / m^2)'},...
{'Surface', 'temperature (^oC)'},...
    {'Total annual', 'melt (mm w.e.)'}};
NameFile = sprintf('%s/_4-SEB',OutputFolder);
PlotWeather_CC(data_surf_yr,vis,VarList,LabelList,NameFile,model_list2);

VarList = {'snowfall','sublimation_mweq','rainfall','SMB_mweq'};
LabelList =  { {'Total annual', 'snowfall (mm w.e.)'},...
    {'Total annual', 'sublimation (mm w.e.)'} ,...
    {'Total annual', 'rainfall (mm w.e.)'},...
    {'Annual surface mass', 'balance (mm w.e.)'} };
if sum(data_surf{4}. runoff)>0.001
    VarList{length(VarList)+1} = 'runoff';
    LabelList{length(LabelList)+1} = 'Total annual \newline runoff (mm w.e.)';
end

NameFile = sprintf('%s/_5-SMB',OutputFolder);
[f,ha] = PlotWeather_CC(data_surf_yr,vis,VarList,LabelList,NameFile,model_list2,1000);
ha(3).Visible = 'off';
cla(ha(3))
ha(4).Position(1) = ha(1).Position(1);
ha(4).Position(3) = 0.8;
ha(4).YAxisLocation='left';
set(f,'CurrentAxes',ha(4))
ha(4).Children(1).String='C';
ha(4).Children(1).Position(2)=620;

sh1 = stairs(datenum(accum_all.year,1,1),accum_all.max*1000,'w','LineStyle','none');
sh2 = stairs(datenum(accum_all.year,1,1),accum_all.min*1000,'w','LineStyle','none');
x = [sh1.XData(1),repelem(sh1.XData(2:end),2)];
y1 = [repelem(sh1.YData(1:end-1),2),sh1.YData(end)];
y2 = [repelem(sh2.YData(1:end-1),2),sh2.YData(end)];
x2 = [x, fliplr(x)];
inBetween = [y1, fliplr(y2)];
h_d =  plot(NaN,NaN,'w');
h_fill =  fill(x2, inBetween, ...
    RGB('light gray'),...
    'LineStyle','none',...
    'FaceAlpha',0.7);
h_med = stairs(datenum(accum_all.year,1,1),accum_all.median*1000,...
    'k','LineWidth',2);
uistack(sh1,'bottom')
uistack(sh2,'bottom')
uistack(h_fill,'bottom')
ylim([80 710])

legendflex([h_fill, h_med],...
    {'min/max enveloppe','median'},...
    'title','Firn core observations',...
    'buffer',[-200 -1],...
    'FontSize',10,'ncol',2);
print(f,NameFile,'-djpeg'); 

%% Comparing melt with racmo and mar
filename = '..\..\Data\RCM\RACMO\Data_AWS_RACMO2.3p2_FGRN055_1957-2017\RACMO_3h_AWS_sites.nc';
time_racmo = ncread(filename,'time')+datenum(1900,1,1);
melt = ncread(filename,'melt');
melt = melt(:,27)*3;
melt(time_racmo<data_surf{1}.time(1)) = [];
time_racmo(time_racmo<data_surf{1}.time(1)) = [];

load('MAR_CC.mat')
data_MAR(data_MAR.time<data_surf{1}.time(1),:) = [];

figure
hold on
plot(data_surf{1}.time,cumsum(data_surf{1}.melt_mweq*1000)-384,'LineWidth',1.5)
plot(data_surf{2}.time,cumsum(data_surf{2}.melt_mweq*1000)-154,'LineWidth',1.5)
plot(data_surf{2}.time,cumsum(data_surf{3}.melt_mweq*1000)-121,'LineWidth',1.5)
plot(data_surf{2}.time,cumsum(data_surf{4}.melt_mweq*1000)-95,'LineWidth',1.5)
plot(time_racmo,cumsum(melt)-175,'LineWidth',1.5)
plot(data_MAR.time,cumsum(data_MAR.ME)-154,'LineWidth',1.5)
plot(data_surf{1}.time,data_surf{1}.time*0,'k')
xlim(data_surf{1}.time([1 end]))
set_monthly_tick(data_surf{1}.time)
ylabel('Cumulated melt (mm w.e.)')
xlabel('Years')
legend('RACMO-adjusted',  'CanESM-adjusted RCP2.6',...
    'CanESM-adjusted RCP4.5', 'CanESM-adjusted RCP8.5',...
    'RACMO2.3p2','MARv3.9.5','Location','northwest',...
    'Interpreter','none')

%% Loading aws data
% addpath(genpath('../AWS_Processing/'))
[data_AWS_GITS] = ImportGCnetData('GITS','yes');
[data_AWS_CEN] = ImportPROMICEData('CEN');

ind = data_AWS_CEN.time<=datenum('31-Dec-2019 21:00:00');
data_AWS_CEN=data_AWS_CEN(ind,:);
            
data_AWS_GITS = TreatAndFilterData(data_AWS_GITS,'ConvertHumidity', 'GITS', OutputFolder, 'off');
data_AWS_CEN = TreatAndFilterData(data_AWS_CEN,'ConvertHumidity', 'CEN', OutputFolder, 'off');
% data_AWS_CEN = ResampleTable(data_AWS_CEN);
data_AWS_CEN. time = data_AWS_CEN.time(1) + (0:length(data_AWS_CEN.time)-1)'/24;
% data_AWS_GITS = ResampleTable(data_AWS_GITS);
data_AWS_GITS. time = data_AWS_GITS.time(1) + (0:length(data_AWS_GITS.time)-1)'/24;

data_AWS_CEN.Albedo = data_AWS_CEN.ShortwaveRadiationUpWm2./data_AWS_CEN.ShortwaveRadiationDownWm2;
data_AWS_CEN.Albedo(data_AWS_CEN.Albedo>1) = NaN;
data_AWS_CEN.Albedo(data_AWS_CEN.Albedo<0.4) = NaN;
data_AWS_GITS.time = data_AWS_GITS.time -1;

data_surf{1}.SRin(2:end) = data_surf{1}.SRin(1:end-1);
data_surf{1}.SRout(2:end) = data_surf{1}.SRout(1:end-1);
data_surf{2}.SRin(1:end-1) = data_surf{2}.SRin(2:end);
data_surf{2}.SRout(1:end-1) = data_surf{2}.SRout(2:end);

%% Reporting available observation data
var_list = {'AirTemperature2C', 'RelativeHumidity2Perc',...
    'WindSpeed2ms','AirPressurehPa',...
    'ShortwaveRadiationDownWm2', 'LongwaveRadiationDownWm2'};
    
for i = 1:length(var_list)
    disp(var_list{i})
    try disp((sum(~isnan(data_AWS_CEN.(var_list{i})))...
        +sum(~isnan(data_AWS_GITS.(var_list{i}))))/24/365)
    catch me
        disp((sum(~isnan(data_AWS_CEN.(var_list{i}))))/24/365)
    end
end
%%  Comparison model vs station data

for ii =1:4
    data_surf_d{ii} = AvgTable(data_surf{ii},'daily',@mean);
end
data_AWS_GITS_d = AvgTable(data_AWS_GITS,'daily',@mean);
data_AWS_CEN_d = AvgTable(data_AWS_CEN,'daily',@mean);

%%
%  close all
VarList = {'theta_2m','RH_2m_wrtw','ws_10m','SRin','LRin'};
VarList2 = {'AirTemperature2C','RelativeHumidity2Perc',...
    'WindSpeed2ms',...
    'ShortwaveRadiationDownWm2',...
    'LongwaveRadiationDownWm2'};
LabelList =  {'Air temperature (^oC)',...
    'Relative humidity (%)',...
    'Wind speed (m/s)',...
    'Downward SW radiation (W/m^2)',...
    'Downward LW radiation (W/m^2)'};
NameFile = sprintf('%s/_2-SurfaceForcing_scatter_after',OutputFolder);
Plot_scatter_CC(data_surf,VarList, data_AWS_GITS,  ...
    data_AWS_CEN,VarList2, vis,LabelList, NameFile,model_list2,'Calibrated')
NameFile = sprintf('%s/2-SurfaceForcing_scatter_after_daily',OutputFolder);
Plot_scatter_CC(data_surf_d,VarList, data_AWS_GITS_d,  ...
    data_AWS_CEN_d,VarList2, vis,LabelList, NameFile,model_list2,'Calibrated daily')

%%  Comparison model before adjustment vs station data
close all

data_RCM{1} =  LoadRCMData('GITS','RACMO');
data_RCM{2} =  LoadRCMData('GITS','CanESM_rcp26');
data_RCM{3} =  LoadRCMData('GITS','CanESM_rcp45');
data_RCM{4} =  LoadRCMData('GITS','CanESM_rcp85');

for i = 1:4

     data_RCM{i} = data_RCM{i}(...
         and(data_RCM{i}.time>=data_surf{i}.time(1),...
         data_RCM{i}.time<=data_surf{i}.time(end)),:);

    data_RCM{i}.WindSpeed2ms=max(0,data_RCM{i}.WindSpeed2ms);     
    time_start = data_surf{i}.time(1);
    time_end = data_surf{i}.time(end);
        % here the RCM file is cropped to the desired period 
    ind_RCM = and(data_RCM{i}.time>=time_start, data_RCM{i}.time<=time_end);
    data_RCM{i} = data_RCM{i}(ind_RCM,:);

    ind_RCM = find(and(data_RCM{i}.time>=time_start, data_RCM{i}.time<=time_end));
    if ~isempty(ind_RCM)
        if time_end>data_RCM{i}.time(end)
            added_rows  = array2table([[data_RCM{i}.time(end)+1/24:1/24:time_end]',...
                NaN( length([data_RCM{i}.time(end)+1/24:1/24:time_end]),size(data_RCM{i},2)-1)],...
                'VariableNames',data_RCM{i}.Properties.VariableNames);
            data_RCM{i} = [data_RCM{i};  added_rows];
        end
    end
     data_RCM_d{i} = AvgTable( data_RCM{i},'daily',@mean);
end

% close all

VarList2 = {'AirTemperature2C','RelativeHumidity2Perc',...
    'WindSpeed2ms',...
    'ShortwaveRadiationDownWm2',...
    'LongwaveRadiationDownWm2'};

LabelList =  {'Air temperature (^oC)',...
    'Relative humidity (%)',...
    'Wind speed (m/s)',...
    'Downward SW radiation (W/m^2)',...
    'Downward LW radiation (W/m^2)'};
NameFile = sprintf('%s/S3-SurfaceForcing_scatter_before',OutputFolder);
Plot_scatter_CC(data_RCM,VarList2, data_AWS_GITS,  ...
    data_AWS_CEN,VarList2, vis,LabelList, NameFile,model_list2,'Raw');
NameFile = sprintf('%s/S3-SurfaceForcing_scatter_before_daily',OutputFolder);
Plot_scatter_CC(data_RCM_d,VarList2, data_AWS_GITS_d,  ...
    data_AWS_CEN_d,VarList2, vis,LabelList, NameFile,model_list2,'Raw daily');

%% Evaluation snowfall before/after adjustment

disp('Adjustement of snowfall rates (%)')
for i = 1:4
    disp(model_list{i})
     disp(nansum(data_surf{i}.snowfall-data_RCM{i}.Snowfallmweq)./...
         nansum(data_RCM{i}.Snowfallmweq)*100);
end

%%  Preparing temperature data
%       Loading RACMO version for depth calculation
disp('Loading racmo run full res')
tic
    % extract run parameters
    load(strcat(path_list{1},'/run_param.mat'));
    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',path_list{1},1);
    finfo = ncinfo(namefile);
%     names={finfo.Variables.Name};
    names = {'time','Year','Day','Hour','SMB_mweq','melt_mweq','H_comp',...
        'runoff','snowfall','rainfall',...
        'sublimation_mweq','Tsurf'};
    
    data_racmo = table();
    for i= 1:length(names)
        eval(sprintf('data_racmo.%s = ncread(''%s'',''%s'');',...
            char(names{i}), namefile,char(names{i})));
    end

    data_racmo.time = ...
        data_racmo.time + datenum(1900,1,1);
data_racmo.time = data_racmo.time(1) + ((1:length(data_racmo.time))'-1)/24;

    % extract subsurface variables
    varname= {'compaction'  'rho' 'slwc' 'snowc' 'snic' 'T_ice'};
    for i =1:length(varname)
        namefile = sprintf('%s/%s_bin_%i.nc',path_list{1},varname{i},1);
        eval(sprintf('%s = ncread(''%s'',''%s'');', varname{i}, namefile,varname{i}));
    end
    fprintf('\nData extracted from nc files.\n');

    % Update BV2017: (not used anymore) Liquid water does not participate to the actual thickness
    % of a layer. Ice does not participate as long as there is enough pore
    % space in the snow to accomodate it.
    pore_space = snowc .* c.rho_water.*( 1./rho - 1/c.rho_ice);
    excess_ice = max(0, snic * c.rho_water / c.rho_ice - pore_space);
    thickness_act = snowc.*(c.rho_water./rho) + excess_ice;
        
    % depth scale from surface
    depth_act_s=zeros(size(thickness_act));
    for j=1:length(thickness_act(1,:))
            depth_act_s(1:end,j)=cumsum(thickness_act(1:end,j));
    end
    data_racmo.H_surf = depth_act_s(end,:)'-depth_act_s(end,1); %+snowbkt*1000/315;

    for i = 1:length(data_racmo.H_surf)-1
        if (data_racmo.H_surf(i+1)-data_racmo.H_surf(i))> c.new_bottom_lay-1
            data_racmo.H_surf(i+1:end) = data_racmo.H_surf(i+1:end) - c.new_bottom_lay*c.rho_water/c.rho_ice;
        end
    end
    
    %depth scale from bottom of the model
    depth_act_m = vertcat(zeros(size(depth_act_s(1,:))), depth_act_s);

    for i = 1:size(depth_act_m,1)
        depth_act_m(i,:) = depth_act_m(i,:) - data_racmo.H_surf' + ...
            data_racmo.H_surf(1);
    end
toc

%       Loading GITS firn temperature 
disp('Loading GITS subsurface temperatures')

c.verbose = 1;

c_all{1}.InputAWSFile = 'C:/Data_save/CC/AWS/GITS_RACMO/data_GITS_combined_hour.txt';
clearvars data_GITS data data_out
[~, ~, ~, ~, ~,...
    ~, ~, ~, ~, ~,~, ...
    ~, ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ~, ~,...
    ~, ~, ~, ~, T_ice_obs, ...
    depth_thermistor, Surface_Height, ~, data_AWS, ~] = ...
    ExtractAWSData(c_all{1});
data_GITS.time = data_AWS.time;
if length(data_AWS.time)==length(data_GITS.time)
    ind = 1:length(data_GITS.time);
else
    ind = and(data_GITS.time<=data_racmo.time(end),...
        data_GITS.time>=data_racmo.time(1));
end
data_GITS.depth = depth_thermistor(ind,:)';
data_GITS.depth(data_GITS.depth==0) = NaN;

data_GITS.T_ice_obs = T_ice_obs(ind,:)' + 273.15;
data_GITS.T_ice_obs(isnan(data_GITS.depth)) = NaN;

data_GITS.T_ice_obs(:,1:(365*2))=NaN;

data_GITS.time = data_GITS.time(ind);
data_GITS.TT_obs= repmat(data_GITS.time',size(data_GITS.depth,1),1);

% Removal of erroneous data from thermistor
 disp ('Removing erroneous data')
 tic
VarName={};
for i =1:32
    VarName{i} = sprintf('IceTemperature%i',i);
end

for i =1:size(data_GITS.T_ice_obs,1)
     data_GITS.T_ice_obs (i, movvar(data_GITS.T_ice_obs(i,:),24*7,'omitnan')>0.05) = NaN;
     data_GITS.T_ice_obs (i, movvar(data_GITS.T_ice_obs(i,:),6,'omitnan')>0.1) = NaN;
%         data_obs.T_ice_obs_out (i, :) = hampel(t_save,24*4,0.1);      
end

data = array2table(data_GITS.T_ice_obs',...
    'VariableNames',{VarName{1:size(data_GITS.T_ice_obs,1)}});
data.time = data_GITS.time;
data_out= SetErrorDatatoNAN(data, 'GITS',vis,OutputFolder,'plot');
data_GITS.T_ice_obs_out = table2array(data_out(:,1:size(data_GITS.T_ice_obs,1)))';

data_GITS.T_ice_obs_sav = data_GITS.T_ice_obs;
data_GITS.T_ice_obs = data_GITS.T_ice_obs_out;
clearvars data_obs.T_ice_obs_sav data_obs.T_ice_obs_out data_out 
toc

% Correcting thermistor location for the GC-Net temp
% the thermistor depth reported in the AWS file cannot be used because of
% missing measured surface height, we therefore need to reprocess it
data_GITS.depth = data_GITS.depth*NaN;
maintenance = ImportMaintenanceFile('GITS');
maintenance(maintenance.date>data_racmo.time(end),:) = [];
maintenance(maintenance.date<data_racmo.time(1),:) = [];
ind_therm_fields = find(contains(maintenance.Properties.VariableNames,'Depth'));
ind_no_new_depth = isnan(nanmean(table2array(maintenance(:,ind_therm_fields)),2));
maintenance(ind_no_new_depth,:) = [];

% for each maintenance entry
for i = 1:size(maintenance,1)           
        % time steps for which the depth will be updated
        ind_time = find(data_GITS.time >= maintenance.date(i));
        data_GITS.depth(:,ind_time) = ...
            repmat(table2array(maintenance(i,ind_therm_fields))',1,length(ind_time))...
            + data_racmo.H_surf(ind_time)' - data_racmo.H_surf(ind_time(1));
end
   
% depth_save = data_obs.depth;

[data_GITS.depth_m, data_GITS.T_ice_obs] = ...
    CalculateThermistorDepth(data_racmo.time,...
    depth_act_s, data_racmo.H_surf,...
    compaction, data_GITS.depth,data_GITS.T_ice_obs_sav, 'GITS');

data_GITS.depth_s=data_GITS.depth_m;
data_GITS.depth_m = repmat(data_racmo.H_surf',10,1)-data_GITS.depth_s;

%       Loading extra thermistor at Camp Century
disp('Loading CEN_THM')
tic

data_THM = struct();
[data_THM.time,~,data_THM.depth, data_THM.T_ice]= LoadCENTHMData();
data_THM.T_ice(1:3,:) = NaN;
data_THM.time_org = data_THM.time;
% cropping/extending the observations to match size of model outputs.
ind_transition = (data_THM.time(1)-data_racmo.time(1))*24;
ind_end = find(data_THM.time(end)==data_racmo.time);
missing_end = length(data_racmo.time) - ind_end;

data_THM.T_ice = [NaN(size(data_THM.T_ice,1),...
    (data_THM.time(1)-data_racmo.time(1))*24),...
    data_THM.T_ice, ...
    NaN(size(data_THM.T_ice,1),missing_end)] ;

data_THM.depth = [NaN(size(data_THM.T_ice,1),...
    (data_THM.time(1)-data_racmo.time(1))*24),...
    data_THM.depth, ...
    NaN(size(data_THM.T_ice,1),missing_end)] ;

data_THM.time = [data_racmo.time(1:ind_transition)',...
    data_THM.time, ...
    data_THM.time(end)+1/24*(1:missing_end)] ;
TT_THM = repmat(data_THM.time,size(data_THM.T_ice,1),1);
clearvars ind_end ind_transition
toc

% Correcting thermistor location for the CEN_THM temp
% addpath(genpath('../../Data/AWS/'))
time_mod = data_racmo.time;

[data_THM.depth_m, data_THM.T_ice_m] = ...
    CalculateThermistorDepth(time_mod,...
    depth_act_s, data_racmo.H_surf,...
    compaction, data_THM.depth, data_THM.T_ice, 'CEN_THM');

data_THM.depth_s=data_THM.depth_m;
data_THM.depth_m = repmat(data_racmo.H_surf',50,1)-data_THM.depth_m;

% Plotting observed temp at GITS & CEN_THM and T_diff
tmp = nansum(~isnan(data_GITS.T_ice_obs),1)>0;
ind1 = find(tmp,1,'first');
ind2 = find(tmp,1,'last');
f = figure('Visible',vis);
subplot(1,2,1)
col1 = PlotTemp(data_GITS.TT_obs(:,ind1:24:ind2),...
    -data_GITS.depth_m(:, ind1:24:ind2),...
    data_GITS.T_ice_obs(:, ind1:24:ind2)-273.15,...
    'PlotTherm', 'yes',...
    'PlotIsoTherm', 'no',...
    'ShowLegend','no',...
    'cmap','jet',...
    'Interp','on',...
    'XLabel','',...
    'YLabel',' ',...
    'CLabel','Firn temperature (^oC)',...
    'Range', -40:1:0,...
    'FlatSurface','no'); 

plot(time_mod(ind1:ind2), -data_racmo.H_surf(ind1:ind2),'LineWidth',2)
axis tight
ylabel('Height above initial surface level (m)','Interpreter','tex')
title('GITS')
subplot(1,2,2)
step  = 72;
ind_time = find(and(data_THM.time>=data_THM.time_org(1),data_THM.time<=data_racmo.time(end)));
col2 = PlotTemp(TT_THM(:,ind_time(1:step:end)),...
    -data_THM.depth_m(:, ind_time(1:step:end)),...
    data_THM.T_ice(:, ind_time(1:step:end)),...
    'PlotTherm', 'yes',...
    'PlotIsoTherm', 'no',...
    'ShowLegend','no',...
    'cmap','jet',...
    'Interp','on',...
    'XLabel','',...
    'YLabel',' ',...
    'CLabel','Firn temperature (^oC)',...
    'Range', -40:1:0,...
    'FlatSurface','no'); 
plot(time_mod(ind_time), -data_racmo.H_surf(ind_time),'LineWidth',2)
plot(time_mod(ind_time), -data_racmo.H_surf(ind_time)+10,'--r','LineWidth',2)
axis tight
title('CEN_THM')
set_monthly_tick(time_mod(ind_time))
datetick('x','mmm-yy','keepticks','keeplimits')
set(gca,'XTickLabelRotation',45)
for i = 2:2:length(col1.YTickLabel)
    col1.YTickLabel(i,:) = '   ';
    col2.YTickLabel(i,:) = '   ';
end
print(f,'./Output/Camp Century/T_firn_obs','-djpeg')
data_THM_save = data_THM;

% 3.6 10 m temp
    disp('Calculating 10 m temperature')
    
% In 1977, firn density was measured to a depth of 100.0 m, and a firn temperature of -24.29°C was measured at 10 m depth (Clausen and Hammer, 1988)
% 5/21/13,B 2-020,77.21666667,61.02166667,1905.8,9.5m,-23,-0.1,-23.1,-23.29,-23.9,0.61

hist_temp = table();
hist_temp.date = ['21-May-2013'; '01-Jul-1977'; '01-Jul-1952'];
hist_temp.time = datenum(hist_temp.date);
hist_temp.T_firn = [-23; -24.29; -23.29];

% 
% GITS obs
ind_nonan = ~isnan(data_GITS.T_ice_obs);
data_GITS.T10 = NaN * data_GITS.time;
for i = 1:length(data_GITS.time)
    if sum(ind_nonan(:,i))>3
        [~, ind_deepest] = max(data_GITS.depth_s(ind_nonan(:,i),i));
        data_GITS.T10(i) = data_GITS.T_ice_obs(ind_deepest,i)-273.15;
    end
end

data_GITS.year = table();
data_GITS.year.time = [7.3094e+05, 7.3094e+05+datenum(1,1,1) ...
    7.3279e+05 7.3279e+05+datenum(1,1,1)];
data_GITS.year.T10 = NaN*data_GITS.year.time;
for i =1:length(data_GITS.year.time)-1
    ind = data_GITS.time==data_GITS.year.time(i);
    data_GITS.year.T10(i) = nanmean_90(data_GITS.T10(ind));
end

% THM obs
[depth_alt, M_alt, ~] = ...
    InterpGridColumnWise(data_THM.depth_s, ...
    data_THM.T_ice, 10, ...
    'lin','VolumeWeighted');
for i = 1:length(M_alt(1,:))
    if sum(~isnan(M_alt(:,i))) == 0
        continue
    end
    M_alt(isnan(M_alt(:,i)),i) = interp1(find(~isnan(M_alt(:,i))),M_alt(~isnan(M_alt(:,i)),i),find(isnan(M_alt(:,i))));
end
data_THM.T10 = M_alt(depth_alt==10);
data_THM.year = table();
ind_y1 = find(~isnan(data_THM.T10),1,'first');
data_THM.year.time = [data_THM.time(ind_y1) data_THM.time(end)];
data_THM.year.T10 = [data_THM.T10(ind_y1), data_THM.T10(end)];

clearvars ind_y1 ind_y2 ind

% simulation
for ii = 1:4
        % temperature at 10 m
    [depth_alt, M_alt, ~] = ...
        InterpGridColumnWise(data_subsurf{ii}.depth_act, ...
        data_subsurf{ii}.T_ice, 10, ...
        'lin','VolumeWeighted');
    data_surf{ii}.T10 = M_alt(depth_alt==10)-273.15;
end

%% 5 plotting firn temperature and 10 m temp
disp('Plotting modelled temperature ')
col = lines(4);
step = 50;
ylimit = 15;
clearvars h h_obs
f = figure('Visible', vis,'OuterPosition',[0 0 23 25]);
set(f,'DefaultAxesFontSize',14)
ha = tight_subplot(8,3,[0.02 0.01],[0.1 0.05],[0.085 0.13]);
for ii =1:4
    set(f,'CurrentAxes',ha(3*(ii-1)+1))
    col_bar = PlotTemp([data_subsurf{ii}.TT(1,1:step:end); data_subsurf{ii}.TT(:,1:step:end)],...
        [-data_surf{ii}.H_surf(1:step:end)'*0;  data_subsurf{ii}.depth_act(:,1:step:end)], ...
        [data_subsurf{ii}.T_ice(1,1:step:end); data_subsurf{ii}.T_ice(:,1:step:end)]-273.15,...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'XLabel','Year',...
        'YLabel','Depth (m)',...
        'CLabel','Modelled firn temperature (^oC)',...
        'cmap','parula',...
        'FlatSurface','no',...
        'Range', -30:2:0);
    plot(data_surf{ii}.time,data_surf{ii}.time*0+10,'--k');

    ylim([0 ylimit])
    col_bar.TickLength=col_bar.TickLength*2;
    set_monthly_tick(data_surf{4}.time)
    if ii==4
        col_bar.Position = [ .9    0.53    0.02    0.42];
    else
        col_bar.Position(1) = 10;
    end
    
    if ii == 1
        set(gca,'XAxisLocation','top')
    else
        set(gca,'XTickLabel','')
    end
    xlim(data_surf{4}.time([1 end]))
    set(gca,'Color',[0.95 0.95 0.95]);
    xlabel('')

    h_now = gca;
    h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
        h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.17, [ABC(ii) ': ' model_list2{ii}]);
    h_text.FontWeight='bold';
    h_text.FontSize=14;
    h_text.Color='w';
end

for i = [2 3 5 6 8 9 11 12 14:21 ]
    ha(i).Visible = 'off';
end
for i = [1 4 7 10 13]
    ha(i).Position(3)= 0.78;
end


col_lines = RGB('gray');
dates =  datenum(['15-Apr-2002'; '23-Apr-2006'; '01-Aug-2017']);

ha(14).Visible = 'off';
ha(14).Position = [0 0 1 1];
set(f,'CurrentAxes',ha(14))
hold on
plot([0.298 0.2 ],[ 0.413 0.37],'LineWidth', 2, 'Color',col_lines);
plot([0.305 0.5 ],[ 0.413 0.37],'LineWidth', 2, 'Color',col_lines);
plot([0.356 0.75 ],[ 0.413 0.37],'LineWidth', 2, 'Color',col_lines);
xlim([0 1])
ylim([0 1])
ha(14).XColor='w';
ha(14).YColor='w';

for i = [1 4 7 10 13]
    set(f,'CurrentAxes',ha(i))
    if i == 13
        hold on
        y=-[15 27];
    else
        y=[0 15];
    end
    h(1) = plot(datenum(dates(1))*[1 1], y,'LineWidth',2,'Color',col_lines);
    h(2) = plot(datenum(dates(2))*[1 1], y,'LineWidth',2,'Color',col_lines);
    h(3) = plot(datenum(dates(3))*[1 1], y,'LineWidth',2,'Color',col_lines);
end

ha(22).Position(4)= 0.27;
ha(23).Position(4)= ha(22).Position(4);
ha(24).Position(4)= ha(22).Position(4);
ha(13).Position(2) = 0.41;

set(f,'CurrentAxes',ha(13))
uistack(ha(13), 'top')
for ii =1:4
    h(ii)=plot(data_surf{ii}.time,data_surf{ii}.T10,'LineWidth',2,'Color',col(ii,:));
    hold on
end

h_obs = scatter(hist_temp.time, hist_temp.T_firn,90,'dk','fill');
scatter(data_GITS.year.time, data_GITS.year.T10,90,'dk','fill')
scatter(data_THM.year.time, data_THM.year.T10,90,'dk','fill')
axis fill; grid on;
% xlabel('Year')
axis tight
ylabel({'10 m firn temperature', '(^oC)'},'interpreter','tex')
xlim([data_surf{ii}.time(1) data_surf{ii}.time(end)])
set_monthly_tick(data_surf{ii}.time,gca);
set(gca,'XTickLabelRotation',0)

legendflex([h, h_obs],[model_list2 'Observation'],'interpreter','none',...
    'anchor',{'n', 'n'},...
    'buffer',[0 27],'nrow',1,'Color','none','box','off')
h_now = gca;
h_text = text(h_now.XLim(2)-(h_now.XLim(2)-h_now.XLim(1))*0.05,...
    h_now.YLim(1)+(h_now.YLim(2)-h_now.YLim(1))*0.2, 'E');
h_text.FontWeight='bold';
h_text.FontSize=14;


set(f,'CurrentAxes',ha(22))
col = lines(4);
hold on
[~, ind_obs] = min(abs(dates(1)-data_GITS.time));
plot(data_GITS.T_ice_obs(~isnan(data_GITS.T_ice_obs(:,ind_obs)),ind_obs)-273.15,...
    -data_GITS.depth_s(~isnan(data_GITS.T_ice_obs(:,ind_obs)),ind_obs),...
    'dk','MarkerFaceColor','k',...
    'LineWidth',2)
for ii = 1:4
    [~, ind_obs] = min(abs(dates(1)-data_surf{ii}.time));
    plot(data_subsurf{ii}.T_ice(:,ind_obs)-273.15,...
        -data_subsurf{ii}.depth_act(:,ind_obs),...
        'Color',col(ii,:),...
    'LineWidth',2)
end
ylim([-20 0])
ylabel('Depth (m)')
box on, grid on;
xlim([-35 1])

set(f,'CurrentAxes',ha(23))
col = lines(4);
hold on
[~, ind_obs] = min(abs(dates(2)-data_GITS.time));
plot(data_GITS.T_ice_obs(~isnan(data_GITS.T_ice_obs(:,ind_obs)),ind_obs)-273.15,...
    -data_GITS.depth_s(~isnan(data_GITS.T_ice_obs(:,ind_obs)),ind_obs),...
    'dk','MarkerFaceColor','k',...
    'LineWidth',2)
for ii = 1:4
    [~, ind_obs] = min(abs(dates(2)-data_surf{ii}.time));
    plot(data_subsurf{ii}.T_ice(:,ind_obs)-273.15,...
        -data_subsurf{ii}.depth_act(:,ind_obs),...
        'Color',col(ii,:),...
    'LineWidth',2)
end
ylim([-20 0])
xlabel('Firn temperature (deg ^oC)','Interpreter','tex')
box on, grid on;
xlim([-35 1])

set(f,'CurrentAxes',ha(24))
col = lines(4);
hold on
[~, ind_obs] = min(abs(dates(3)-data_THM.time));
plot(data_THM.T_ice(~isnan(data_THM.T_ice(:,ind_obs)),ind_obs),...
    -data_THM.depth_s(~isnan(data_THM.T_ice(:,ind_obs)),ind_obs),...
    'dk','MarkerFaceColor','k',...
    'LineWidth',2)
for ii = 1:4
    [~, ind_obs] = min(abs(dates(3)-data_surf{ii}.time));
    plot(data_subsurf{ii}.T_ice(:,ind_obs)-273.15,...
        -data_subsurf{ii}.depth_act(:,ind_obs),...
        'Color',col(ii,:),...
    'LineWidth',2)
end
box on, grid on;
ylim([-20 0])
xlim([-35 1])

count = 5;
for i = [22 23 24]
    set(f,'CurrentAxes',ha(i))
    count = count+1;
    h_now = gca;
    h_text = text(h_now.XLim(2)-(h_now.XLim(2)-h_now.XLim(1))*0.1,...
        h_now.YLim(1)+(h_now.YLim(2)-h_now.YLim(1))*0.1, ABC(count));
    h_text.FontWeight='bold';
    h_text.FontSize=14;
    
    if i>22
        set(ha(i),'YTickLabel','')
    end
end


print(f, sprintf('%s/_7-Temp_mod',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% temperature table

temp_obs=table();
temp_obs.time = [hist_temp.time; data_GITS.year.time'; data_THM.year.time'];
temp_obs.T10 = [hist_temp.T_firn; data_GITS.year.T10'; data_THM.year.T10'];
[~, ind] = sort(temp_obs.time);
temp_obs=temp_obs(ind,:);
fprintf("Model; Trend 1966-2020 (^oC a-1); Trend 2021-2100 (^oC a-1);  ME;  RMSE; \n")
trend_10m = [];
for ii =1:4
    lm1 = fitlm(data_surf{ii}.time(data_surf{ii}.time<datenum(2020,1,1)),...
        data_surf{ii}.T10(data_surf{ii}.time<datenum(2020,1,1)));
    if ii==1
        lm2=[];
        lm2.Coefficients.Estimate=[NaN NaN];
    else
        lm2 = fitlm(data_surf{ii}.time(data_surf{ii}.time>=datenum(2020,1,1)),...
            data_surf{ii}.T10(data_surf{ii}.time>=datenum(2020,1,1)));
    end
    
    [LIA] = find(ismember(data_surf{ii}.time,temp_obs.time));
    LOCB= find(ismember(temp_obs.time,data_surf{ii}.time));
    ME = nanmean(data_surf{ii}.T10(LIA)-temp_obs.T10(LOCB));
    RMSE = sqrt(nanmean((temp_obs.T10(LOCB) - data_surf{ii}.T10(LIA)).^2));
    
    fprintf("%s; %0.2f; %0.2f;  %0.2f;  %0.2f; \n",model_list2{ii},...
        lm1.Coefficients.Estimate(2)*365*10,...
        lm2.Coefficients.Estimate(2)*365*10,...
        ME, RMSE)
end
% h_obs = scatter(hist_temp.time, hist_temp.T_firn,90,'dk','fill');
% scatter(data_GITS.year.time, data_GITS.year.T10,90,'dk','fill')
% scatter(data_THM.year.time, data_THM.year.T10,90,'dk','fill')

%% 6 - diff plot
data_THM = data_THM_save;

for ii = 1:4
    [~, ind_end] = min(abs(data_THM.time-data_surf{ii}.time(end)));
    ind_transition = (data_THM.time(1)-data_surf{ii}.time(1))*24;

    data_THM.T_ice = [NaN(size(data_THM.T_ice,1),...
        (data_THM.time(1)-data_surf{ii}.time(1))*24),...
        data_THM.T_ice(:,1:ind_end)] ;

    data_THM.depth = [NaN(size(data_THM.T_ice,1),...
        (data_THM.time(1)-data_surf{ii}.time(1))*24),...
        data_THM.depth(:,1:ind_end)] ;
    
    data_THM.depth_s = [NaN(size(data_THM.T_ice,1),...
        (data_THM.time(1)-data_surf{ii}.time(1))*24),...
        data_THM.depth_s(:,1:ind_end)] ;
     data_THM.depth_m = [NaN(size(data_THM.T_ice,1),...
        (data_THM.time(1)-data_surf{ii}.time(1))*24),...
        data_THM.depth_m(:,1:ind_end)] ;
     data_THM.T_ice_m = [NaN(size(data_THM.T_ice,1),...
        (data_THM.time(1)-data_surf{ii}.time(1))*24),...
        data_THM.T_ice_m(:,1:ind_end)] ;

    data_THM.time = [data_surf{ii}.time(1:ind_transition)',...
        data_THM.time(1:ind_end)] ;
    TT_THM = repmat(data_THM.time,size(data_THM.T_ice,1),1);
    clearvars ind_end ind_transition

    f = figure;
    subplot(2,1,1)
    [T_mod_GITS, T_diff_GITS, ME, RMSE, R2]=T_diff_plotting(data_GITS.TT_obs, data_GITS.depth_s,...
        data_GITS.T_ice_obs-273.15,  data_subsurf{ii}.depth_act, data_surf{ii}.Tsurf, T_ice-273.15, ...
        data_surf{ii}.H_surf*0, model_list2{ii}, 'GITS', -10:1:10,OutputFolder, vis);
    fprintf('%s; GITS; %0.2f; %0.2f; %0.2f\n', model_list2{ii}, ME, RMSE, R2);
    xlabel('')
        h_now = gca;
    h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
        h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.15, 'A');
    h_text.FontWeight='bold';
    h_text.FontSize=16;
    subplot(2,1,2)
    [T_mod_THM, T_diff_THM, ME, RMSE, R2]=T_diff_plotting(TT_THM, data_THM.depth_s,...
        data_THM.T_ice_m,  data_subsurf{ii}.depth_act, data_surf{ii}.Tsurf, T_ice-273.15, ...
        data_surf{ii}.H_surf*0, model_list2{ii}, 'CEN-THM', -10:1:10,OutputFolder, vis);
    fprintf('%s; THM; %0.2f; %0.2f; %0.2f\n', model_list2{ii}, ME, RMSE, R2);
    ylim([0 20])
        h_now = gca;
    h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
        h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.15, 'B');
    h_text.FontWeight='bold';
    h_text.FontSize=16;
    
     print(f, sprintf('%s/6-T_diff_%s',OutputFolder, model_list{ii}), '-djpeg')
end

%%  7 - Density
 disp('Plotting modelled density ')
step = 24;
ylimit = 40;

f = figure('Visible', vis,'OuterPosition',[0 0 23 25]);
num_line = 8;
ha = tight_subplot(num_line,5,[0.02 0.04],[0.06],[0.085 0.07]);
for i = 6:(num_line*5-5)
    ha(i).Visible = 'off';
end
ha(2).Position = ha(6).Position;
ha(3).Position = ha(11).Position;
ha(4).Position = ha(16).Position;
for i =1:4
    ha(i).Position(3) = 0.87;
end

count = 0;
for i =(num_line*5-4):(num_line*5)
    ha(i).Position(1) = ha(i).Position(1)-0.15*count/5;
    ha(i).Position(2) = 0.08;
    ha(i).Position(4) = 0.41;
    count = count+1;
end
ha(5).Position = [0 0 1 1];
set(f,'CurrentAxes',ha(5))
hold on
col_lines = RGB('dark gray');
plot([0.085 0.16 ],[ 0.51 0.49],'LineWidth', 2, 'Color',col_lines);
plot([0.26 0.3 ],[ 0.51 0.49],'LineWidth', 2, 'Color',col_lines);
plot([0.338 0.46 ],[ 0.51 0.49],'LineWidth', 2, 'Color',col_lines);
plot([0.375 0.6 ],[ 0.51 0.49],'LineWidth', 2, 'Color',col_lines);
plot([0.86 0.75 ],[ 0.51 0.49],'LineWidth', 2, 'Color',col_lines);
xlim([0 1])
ylim([0 1])
ha(5).XColor='w';
ha(5).YColor='w';
if sum(data_subsurf{4}.rho_all(data_subsurf{4}.depth_act<30)>830)>10
    cmap_dens = 'parula_black';
else
    cmap_dens= 'parula';
end

for ii =1:4
    set(f,'CurrentAxes',ha(ii))
    uistack(ha(ii),'top')
    colbar = PlotTemp([data_subsurf{ii}.TT(1,1:step:end); data_subsurf{ii}.TT(:,1:step:end)],...
        [0*data_surf{ii}.H_surf(1:step:end)';  data_subsurf{ii}.depth_act(:,1:step:end)], ...
        [data_subsurf{ii}.rho_all(1,1:step:end); data_subsurf{ii}.rho_all(:,1:step:end)],...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'XLabel','Year',...
        'YLabel','Depth (m)',...
        'CLabel','Modelled firn density (kg m^{-3})',...
        'cmap',cmap_dens,...
        'Interp','yes',...
        'FlatSurface','no',...
        'Range', 310:20:850);
    tmp = [1967 1996 2010 2017 2100];
    for k = 1:5
        if ii == 1 && k==5
            continue
        end
    plot(datenum(tmp(k)*[1 1], 1, 1), [0 ylimit],'Color',col_lines,...
        'LineWidth',2)
    end
    ylim([0 ylimit])
    set_monthly_tick(data_surf{4}.time)
    if ii==1
        set(gca,'XAxisLocation','top'); xlabel('Year');
        colbar.Position = [ .9    0.55    0.02    0.4000];
    else
        colbar.Position(1) = 10;
        set(gca,'XTickLabel',''); xlabel('');
    end
    temp = colbar.YTickLabel;
    for i = 1:length(temp)
        if i/2==floor(i/2)
            temp(i,:)=' ';
        end
    end
    colbar.YTickLabel = temp;
    xlim(data_surf{4}.time([1 end]))
    set(gca,'Color',[0.95 0.95 0.95]);
    h_tit = title([ABC(ii) ': '...
        model_list2{ii}]);
    h_tit.Parent = ha(ii); h_tit.Units='normalized';
    if ii==1
        h_tit.Position(1:2) = [0.075 0.013];
    else
        h_tit.Position(1:2) = [0.13 0.013];
    end
end

% Validation with cores
disp('Plotting comparison with density profiles from cores')
ii = 1;
NameStation = 'GITS';
i_core = FindCore(Core,'NearestCodeLocation','Camp Century');  

% ordering cores in chronological order
dates = zeros(size(i_core));
for i = 1:length(i_core)
    dates(i) = datenum(Core{i_core(i)}.Info.DateCored);
end
[~, i_ordered] = sort(dates);
i_core = i_core(i_ordered);

% removing the cores that were drilled more than 5 years before the
% beginning of the model run
i_remove = [];
for i = 1:length(i_core)
    DV = datevec(data_surf{ii}.time(1));
    if Core{i_core(i)}.Info.DateCored.Year < DV(:,1) - 1
        i_remove = [i_remove, i];
    end
end
i_core(i_remove) = [];

count = num_line*5-5;

ylim_core = 90;
count = count+1;
    set(f,'CurrentAxes',ha(count))
    time_core = datenum(1966,1,1);
    hold on
    [~, ind_time] = min(abs(data_surf{1}.time - time_core));
    stairs(data_subsurf{1}.rho_all(:,ind_time),...
        data_subsurf{1}.depth_act(:,ind_time),...
        'LineWidth',2,'Color','k');
    set(gca,'YDir','reverse'); ylim([0 80]); xlim([300 900]); grid on; box on
     ylabel('Depth (m)')
    h_tit = title([ABC(count-num_line*5+9) ': 1966']);
    h_tit.Parent = ha(count); h_tit.Units='normalized';
    h_tit.Position(1:2) = [0.3 0.01];
    
for ii = i_core
    count = count+1;
    set(f,'CurrentAxes',ha(count))
    time_core = datenum(Core{ii}.Info.DateCored);

col = lines(4);
    hold on
    for kk = 1:4
        [~, ind_time] = min(abs(data_surf{kk}.time - time_core));
        stairs(data_subsurf{kk}.rho_all(:,ind_time),...
            data_subsurf{kk}.depth_act(:,ind_time),...
            'LineWidth',2,'Color',col(kk,:));
        h_tit = title([ABC(count-num_line*5+9) ': '...
            num2str(Core{ii}.Info.DateCored.Year)]);
        h_tit.Parent = ha(count); h_tit.Units='normalized';
        h_tit.Position(1:2) = [0.2 0.01];
    end
    plot(Core{ii}.Data.Density,Core{ii}.Data.Depth/100,'k','LineWidth',2)

    set(gca,'YDir','reverse','YTickLabel',' ');
    ylim([0 80]); xlim([300 900]); grid on; box on

    h_tit = title([ABC(count-num_line*5+9) ': '...
        num2str(Core{ii}.Info.DateCored.Year)]);
    h_tit.Parent = ha(count); h_tit.Units='normalized';
    h_tit.Position(1:2) = [0.3 0.01];
    if count == 38
        xlabel('Firn density (kg m^{-3})','Interpreter','tex')
    end
end
legendflex([ model_list2, 'Observed'],...
    'ref',gcf,...
    'anchor',{'e','e'},...
    'buffer',[-2 -150],...
    'fontsize',11,...
    'ncol',1,...
    'interpreter','none',...
    'box','off');

 count = count+1;
    set(f,'CurrentAxes',ha(count))
    time_core = datenum(2100,1,1);
    hold on
    for kk = 2:4
        [~, ind_time] = min(abs(data_surf{kk}.time - time_core));
        stairs(data_subsurf{kk}.rho_all(:,ind_time),...
            data_subsurf{kk}.depth_act(:,ind_time),...
            'LineWidth',2,'Color',col(kk,:));
    end

    set(gca,'YDir','reverse','YTickLabel',' ');
    ylim([0 80]); xlim([300 900]); grid on; box on

    h_tit = title([ABC(count-num_line*5+9) ': 2100']);
    h_tit.Parent = ha(count); h_tit.Units='normalized';
    h_tit.Position(1:2) = [0.3 0.01];

    print(f, sprintf('%s/_8-rho_mod',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end
%%
f = figure('OuterPosition',[0 0 8 15]);
col_parula = hsv(4);% brewermap(4,'Paired');
% col_parula([2 4 6],:) = [];
hold on
name={};
count = 1;
h=[];
h(1)=plot(data_subsurf{1}.rho_all(:,1),...
        data_subsurf{1}.depth_act(:,1),...
        'LineWidth',2,'Color',col_parula(1,:));

for ii = (i_core)
    count = count+1;
        h(count) = plot(Core{ii}.Data.Density,Core{ii}.Data.Depth/100,...
            'LineWidth',2,'Color',col_parula(count,:));
end

name = {'1966','1996','2010','2017'};
uistack(h(2),'top')
uistack(h(3),'top')
uistack(h(1),'top')
% plot(data_subsurf{1}.rho_all(:,1)+40,...
%         data_subsurf{1}.depth_act(:,1),...
%         'LineWidth',2,'Color',col_parula(1,:));
% plot(data_subsurf{1}.rho_all(:,1)-40,...
%         data_subsurf{1}.depth_act(:,1),...
%         'LineWidth',2,'Color',col_parula(1,:));
legend(h,name,'interpreter','tex','location','SouthWest','FontSize',10)
set(gca,'YDir','reverse')
axis tight; box on;
xlabel('Firn density (kg m^{-3})','Interpreter','tex')
ylabel('Depth (m)')
    print(f, sprintf('%s/7-Core',OutputFolder), '-djpeg')

%%  8 - Compaction 
opts = delimitedTextImportOptions("NumVariables", 19);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["TRANSMISSION_TIMESTAMP", "YEAR", "MONTH", "DAY", "BATTERY_MIN_V", "BATTERY_MAX_V", "MEAN_PANEL_TEMP_C", "INST_1_RATIO", "INST_1_CORRECTION_RATIO", "INST_2_RATIO", "INST_2_CORRECTION_RATIO", "INST_3_RATIO", "INST_3_CORRECTION_RATIO", "INST_1_LENGTH_UNCORRECTED_M", "INST_2_LENGTH_UNCORRECTED_M", "INST_3_LENGTH_UNCORRECTED_M", "INST_1_LENGTH_CORRECTED_M", "INST_2_LENGTH_CORRECTED_M", "INST_3_LENGTH_CORRECTED_M"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
data_comp = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Data\Camp Century\Compaction\FIRNCO_DATA_2020.01.20.csv", opts);
clear opts
data_comp.time = datenum(data_comp.YEAR,data_comp.MONTH,data_comp.DAY);
data_comp(:,5:16) = [];

% The three useful columns are the last three: "INST_1_LENGTH_CORRECTED_M" 
% (and INST_2, INST_3, respectively). This shows the length of the 2-m 
% potentiometer cable over time. 

% CURRENT_BOREHOLE_LEN = INIT_BOREHOLE_LEN - INIT_CABLE_LEN + CURRENT_CABLE_LEN
ini_top = [1.4 0 0 0];
ini_bot = [62.30 20 5 4.9];
%    4.9 7.7 NaN];
% the two borehole depths were approximately 4.9 +/- 0.2 m and 7.7 +/- 0.2 m.

time_maintenance = [datenum('01-Aug-2017') datenum('16-May-2019')];

[~, ind_stop] = min(abs(data_surf{1}.time - data_comp.time(end)));
data_comp.B_len = NaN(length(data_comp.time), length(ini_top));
for ii = 1:4
    for j=1:4
        B{j,ii} = [data_surf{ii}.H_surf * NaN,  data_surf{ii}.H_surf * NaN];
    end
end
settling_time = 30;
for j = 1:4
    if j<4
        ind_start = find(data_comp.time==time_maintenance(1));
        ind_end = find(data_comp.time==time_maintenance(2)-3);
        data_comp.B_len(ind_start:ind_end,j) = (ini_bot(j)-ini_top(j)) ...
            - data_comp.(['INST_' num2str(j) '_LENGTH_CORRECTED_M'])(ind_start) + ...
            data_comp.(['INST_' num2str(j) '_LENGTH_CORRECTED_M'])(ind_start:ind_end);
        data_comp.B_len(ind_start:(ind_start+settling_time),j) = NaN;

    else
        [~, ind_start] = min(abs(data_comp.time-time_maintenance(2)));
        ind_end = length(data_comp.time);
        data_comp.B_len(ind_start:ind_end,j) = (ini_bot(j)-ini_top(j)) ...
            - data_comp.(['INST_1_LENGTH_CORRECTED_M'])(ind_start) + ...
            data_comp.(['INST_1_LENGTH_CORRECTED_M'])(ind_start:ind_end);
        data_comp.B_len(ind_start:(ind_start+settling_time),j) = NaN;
    end
    
    ini_top(j)=ini_bot(j)-data_comp.B_len(ind_start+settling_time+1,j);

end

step = 12;
time_maintenance = time_maintenance+settling_time;

for ii = 1:4
        
        disp(ii)
        for j = 1:4
            if j<4
                [~, ind_start_mod] = min(abs(data_surf{ii}.time - time_maintenance(1)));
                [~, ind_stop_mod] = min(abs(data_surf{ii}.time - time_maintenance(2)));
            else
                [~, ind_start_mod] = min(abs(data_surf{ii}.time - time_maintenance(2)));
                [~, ind_stop_mod] = min(abs(data_surf{ii}.time - data_comp.time(end)));
            end
            if abs(ind_start_mod-ind_stop_mod)<3
                continue
            end
            [tmp] = track_horizon(ini_top(j),ind_start_mod,step,...
                                        data_subsurf{ii}.depth_act,...
                                        data_subsurf{ii}.compaction, ...
                                        data_surf{ii}.H_surf,...
                                        ind_stop_mod); 
            B{j,ii}(ind_start_mod:end,1) =tmp(ind_start_mod:end);
            
            tmp = track_horizon(ini_bot(j),ind_start_mod,step,...
                                        data_subsurf{ii}.depth_act,...
                                        data_subsurf{ii}.compaction, ...
                                        data_surf{ii}.H_surf,...
                                        ind_stop_mod);
            B{j,ii}(ind_start_mod:end,2) =tmp(ind_start_mod:end);
            
%             B{j,ii}(ind_start_mod:(ind_start_mod+()*24),:) =NaN;
        end
    end
% end

%% Plotting compaction rates and boreholes
 f = figure('OuterPosition',[0 0 24 13],'DefaultAxesFontSize',13);
ha= tight_subplot(2,4,0.05,[0.12 ],[0.09 0.02]);                      

% Comparing borehole length
col = lines(4);
for j = 1:4
    set(f,'CurrentAxes',ha(j))
    hold on
    plot(data_comp.time,data_comp.B_len(:,j),'k','LineWidth',2)
    for ii = 1:4
        [~, ind_end] = min(abs(data_surf{ii}.time - data_comp.time(end)));

        plot(data_surf{ii}.time(1:ind_end),...
            B{j,ii}(1:ind_end,2)-B{j,ii}(1:ind_end,1),...
            'Color',col(ii,:),'LineWidth',2)
    end
    set_monthly_tick(data_comp.time)
    if ismember(j,[1 2 3 5 6 7])
    xlim([time_maintenance(1) datenum(2019,5,1)])
    else
    xlim([time_maintenance(2)+2 data_comp.time(end)])
    end
        
    if j==1
        ylabel('Borehole\newline length (m)','Interpreter','tex')
    end
    h_now = gca;
    h_text = text(h_now.XLim(2)-(h_now.XLim(2)-h_now.XLim(1))*0.2,...
        h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.12, [ABC(j)]);
    h_text.FontWeight='bold';
    h_text.FontSize=16;
end

% Comparing Compaction rates
col = lines(4);
    for j = 1:4
        set(f,'CurrentAxes',ha(4 + j))
        hold on
        fprintf('Instrument %i; ME; RMSE; R2\n',j)
        for ii = 1:4
            tmp=-smooth(gradient(data_comp.B_len(:,j)),7*4)*1000;
            tmp(isnan(data_comp.B_len(:,j)))=NaN;
            B_mod = table();
            B_mod.time = data_surf{ii}.time;
            B_mod.len = B{j,ii}(:,2)-B{j,ii}(:,1);
            B_mod_daily = AvgTable(B_mod,'daily',@mean);
            [~, ind_uni, ~] = unique(B_mod_daily.time);
            B_mod_daily = B_mod_daily(ind_uni,:);
            [~, ind_end] = min(abs(B_mod_daily.time - data_comp.time(end)));
            ind_start = find(B_mod_daily.time==data_comp.time(1));
            tmp1= -(gradient(B_mod_daily.len(ind_start:ind_end)))*1000;
            time_tmp1 = B_mod_daily.time(ind_start:ind_end);
            plot(time_tmp1,...
                tmp1,...
                'Color',col(ii,:),'LineWidth',2)

            tmp = tmp(ismember(data_comp.time,time_tmp1));
            tmp1 = tmp1(ismember(time_tmp1,data_comp.time));
            ME = nanmean(tmp1-tmp);
            RMSE = sqrt(nanmean((tmp-tmp1).^2));
            lm = fitlm(tmp,tmp1);
            fprintf('%s; %0.2f (%0.1f%%); %0.2f (%0.1f%%); %0.2f \n',...
                model_list2{ii},ME,ME/nanmean(tmp(~isnan(tmp1)))*100,...
                RMSE,RMSE/nanmean(tmp(~isnan(tmp1)))*100,...
                lm.Rsquared.Ordinary)

        end
        plot(data_comp.time,tmp,'k','LineWidth',2)
        
        set_monthly_tick(data_comp.time)
        if ismember(j,[1 2 3 5 6 7])
        xlim([time_maintenance(1) datenum(2019,5,1)])
        else
        xlim([time_maintenance(2)+2 data_comp.time(end)])
        end
        
        if j==1
            ylabel('Compaction rate\newline       (mm/day)','Interpreter','tex')
        end
        h_now = gca;
        h_text = text(h_now.XLim(2)-(h_now.XLim(2)-h_now.XLim(1))*0.2,...
            h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.12, [ABC(4 +j)]);
        h_text.FontWeight='bold';
        h_text.FontSize=16;
    end
    
ha(2).XAxisLocation='top';

h_leg = legendflex([model_list2, 'Observation'],...
    'ref',gcf,...
    'anchor',{'n','n'},...
    'buffer',[0 0],...
    'Interpreter','none',...
    'nrow',2,...
    'box','off');
for k = 1:4
    ha(k).XTickLabel = ' ';
end
for k = 1:8
    xti = ha(k).XAxis.MinorTickValues;

    ha(k).LineWidth = 1.5;
    ha(k).TickLength = ha(k).TickLength*2;
    ha(k).XMinorTick='on';
    ha(k).XAxis.MinorTickValues=xti;
    ha(k).TickLength = [0.05 0.05];
    if k>4
    xlabel(ha(k),'Year');
    end
end

print(f, sprintf('%s/_9-compaction_val',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% Calculating   Tracking horizons and Annual layers from cores
clearvars h
disp('Tracking 1966 and 2012 horizon')
for ii =1:4
    disp(ii)
    [~, ind_start] = min(abs(data_surf{ii}.time - datenum(1966,6,1)));
    [data_surf{ii}.depth_66] = track_horizon(0,ind_start,6,...
                                data_subsurf{ii}.depth_act,...
                                data_subsurf{ii}.compaction, ...
                                data_surf{ii}.H_surf);
    [~, ind_start] = min(abs(data_surf{ii}.time - datenum(2012,6,1)));
    [data_surf{ii}.depth_12] = track_horizon(0,ind_start,6,...
                                data_subsurf{ii}.depth_act,...
                                data_subsurf{ii}.compaction, ...
                                data_surf{ii}.H_surf);
    [~, ind_start] = min(abs(data_surf{ii}.time - datenum(2017,6,1)));
    [data_surf{ii}.depth_17] = track_horizon(32,ind_start,6,...
                                data_subsurf{ii}.depth_act,...
                                data_subsurf{ii}.compaction, ...
                                data_surf{ii}.H_surf);
end

% Annual layers from cores
clearvars h
disp('Tracking horizons from JPs files')

% Camp Century 1977 cores 1 and 2 and Camp Century 2010.
% You will find three columns in each file:
% Depth of bottom of annual layer, calendar year of layer and mean delta
% O18 of annual layer (no need of depth of top of layer as cores are continuous).

opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["depth", "year", "mean_dO18"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
CC77_1 = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Data\Camp Century\KU\CC771_annual.txt", opts);
clear opts

opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["depth", "year", "mean_dO18"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
CC77_2 = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Data\Camp Century\KU\CC772_annual.txt", opts);
clear opts

opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["depth", "year", "mean_dO18"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
CC10 = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Data\Camp Century\KU\CC10_annual.txt", opts);
clear opts

years_all =sort( unique([CC77_1.year; CC77_2.year; CC10.year]));
years_all(datenum(years_all,1,1)<data_surf{ii}.time(1)) = [];
depth_annual_layers = cell(1,4);

for ii =1:4
    disp(ii)
    depth_annual_layers{ii} = NaN(length(years_all), length(data_surf{ii}.time));
    for j = 1:length(years_all)
        fprintf('%i / %i\n',j,length(years_all));
        [~, ind_start] = min(abs(data_surf{ii}.time - datenum(years_all(j)-1,9,1)));
        depth_annual_layers{ii}(j,:) = track_horizon(0,ind_start,48,...
                                    data_subsurf{ii}.depth_act,...
                                    data_subsurf{ii}.compaction, ...
                                    data_surf{ii}.H_surf);
    end
end

%% 9 - plotting Tracking horizons and Annual layers from cores
col = lines(4);
    ind_77_in_mod = ismember(CC77_1.year, years_all);
    ind_mod_in_77 = ismember(years_all, CC77_1.year);
    ind_10_in_mod = ismember(CC10.year, years_all);
    ind_mod_in_10 = ismember(years_all, CC10.year);
[~, ind_1977] =  min(abs(data_surf{1}.time - datenum(1977,7,1)));
[~, ind_2010] =  min(abs(data_surf{1}.time - datenum(2010,7,1)));
axis tight


f = figure('Visible',vis,'OuterPosition',[0 0 20 13]);
ha = tight_subplot(1,2,0.1,[0.15 0.1],0.1);
% ha(2).Visible = 'off';
% ha(6).Visible = 'off';
% ha(1).Position(3) = 0.8;
% ha(5).Position(3) = 0.8;

% set(f,'CurrentAxes',ha(1))
% hold on
% for ii =1:1
% %     x = data_surf{ii}.time;
% %     y1 = depth_annual_layers{ii}(find(ind_mod_in_77,1,'first'),:)';
% %     y2 = depth_annual_layers{ii}(find(ind_mod_in_77,1,'last'),:)';
% %     x2 = [x, fliplr(x)];
% %     ind1 = min(find(~isnan(y1),1,'first'), find(~isnan(y2),1,'first'))+10;
% %     ind2 = min(find(~isnan(y1),1,'last'), find(~isnan(y2),1,'last'));
% %     y1=y1(ind1:ind2);
% %     y2=y2(ind1:ind2);
% %     y1(isnan(y1))=0;
% %     y2(isnan(y2))=0;
% %     inBetween = [y1, fliplr(y2)];
% %     x2=x2(ind1:ind2,:);
% %     fill(x2, inBetween, RGB('light light green'));
% 
% plot(data_surf{ii}.time,...
%         depth_annual_layers{ii}([find(ind_mod_in_77,1,'first') ...
%             find(ind_mod_in_77,1,'last')],:),...
%             'Color',RGB('light green'));
%         
% plot(data_surf{ii}.time,...
%         depth_annual_layers{ii}([find(ind_mod_in_10,1,'first') ...
%             find(ind_mod_in_10,1,'last')],:),...
%             'Color',RGB('light blue'));
% end
% plot(data_surf{ii}.time(ind_1977)*ones(size(CC77_1.depth(ind_77_in_mod))),...
%     CC77_1.depth(ind_77_in_mod),'ok','MarkerFaceColor','k')
% plot(data_surf{ii}.time(ind_2010)*ones(size(CC10.depth(ind_10_in_mod))),...
%     CC10.depth(ind_10_in_mod),'ok','MarkerFaceColor','k')
% 
% set(gca,'YDir','reverse','XAxisLocation','top')
% xlabel('Year'); ylabel('Depth (m)');
% set_monthly_tick(data_surf{ii}.time);
% xlim(data_surf{ii}.time([1 end]))
    
set(f,'CurrentAxes',ha(1))
hold on
h = [];
for ii =4:-1:1
    h(ii) = plot(CC77_1.depth(ind_77_in_mod),...
        flipud(depth_annual_layers{ii}(ind_mod_in_77,ind_1977)),'o',...
        'Color',col(ii,:),'MarkerFaceColor',col(ii,:));
    h(ii) = plot(CC77_2.depth(ind_77_in_mod),...
        flipud(depth_annual_layers{ii}(ind_mod_in_77,ind_1977)),'x',...
        'Color',col(ii,:),'MarkerFaceColor',col(ii,:));
end
    b(1) = plot(NaN,NaN,'xk', 'MarkerFaceColor','k');
    b(2) = plot(NaN,NaN,'ok', 'MarkerFaceColor','k');
    text(CC77_1.depth(ind_77_in_mod),...
        flipud(depth_annual_layers{ii}(ind_mod_in_77,ind_1977))+2,...
        num2str(CC77_1.year(ind_77_in_mod)),'Rotation',90);
    plot([0 40],[0 40],'k')
    axis tight; box on; grid on
    xlim([3 10])
    ylim([3 14])
    xlabel('Observed depth (m)')
    ylabel('Simulated depth (m)')
title('Annual layers in the 1977 cores','FontSize',13)
% legend(h,model_list,'location','eastoutside','interpreter','none')
h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'A');
h_text.FontWeight='bold';
h_text.FontSize=16;

set(f,'CurrentAxes',ha(2))
hold on
h = [];
for ii =4:-1:1
    ind1 = find(ismember(CC10.year, years_all));
    ind2 = find(ismember(years_all, CC10.year));
    h(ii) = plot(CC10.depth(ind1),...
        flipud(depth_annual_layers{ii}(ind2,ind_2010)),'o',...
        'Color',col(ii,:),'MarkerFaceColor',col(ii,:));
end
    text(CC10.depth(ind1(1:2:end)),...
        flipud(depth_annual_layers{ii}(ind2((1:2:end)),ind_2010))+2,...
        num2str(CC10.year(ind1((1:2:end)))),'Rotation',90);
    plot([0 40],[0 40],'k')
    axis tight; box on; grid on
    xlim([2 30])
    ylim([2 45])
title('Annual layers in the 2010 core','FontSize',13)
% legend(h,model_list2,'location','eastoutside','interpreter','none')
    xlabel('Observed depth (m)')
    ylabel('Simulated depth (m)')
%     ha(4).YAxisLocation = 'right';
%     ha(3).XMinorTick = 'on';
%     ha(3).XMinorTick = 'on';

    h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'B');
h_text.FontWeight='bold';
h_text.FontSize=16;
print(f, sprintf('%s/9-tracking annual layers',OutputFolder), '-djpeg') 

%% 10- Tracking debris horizon and water percolation
    f = figure('Visible',vis,'OuterPosition',[0 0 28 13]);
h = [];
g=[];
ha= tight_subplot(1,2,0.08,[0.3 0.18] ,0.12 );                      
% ha(2).Position(4) = 0.4;
col = lines(4);
set(f,'CurrentAxes',ha(1))
hold on
for ii =2:4
    h(ii-1) = plot(data_surf{ii}.time,data_surf{ii}.depth_17,...
        'Color',col(ii,:),'LineWidth',2);
%        plot(data_surf{ii}.time,data_surf{ii}.depth_66,...
%         'Color',col(ii,:),'LineWidth',2);
%        plot(data_surf{ii}.time,data_surf{ii}.depth_12,...
%         'Color',col(ii,:),'LineWidth',2);
       
    depth_1 = data_subsurf{ii}.depth_act;
    depth_1(data_subsurf{ii}.slwc<1e-6) = 0;
    perc_depth = max(depth_1,[],1);
    
%     g(ii-1) = plot(data_surf{ii}.time,...
%         perc_depth,'-',...
%         'Color',col(ii,:),'LineWidth',2);
    set(f,'CurrentAxes',ha(2))
    hold on
    g(ii-1) = plot(data_surf{ii}.time,...
        perc_depth,'-',...
        'Color',col(ii,:),'LineWidth',1);
    disp(model_list{ii})
    disp(max(perc_depth))
    disp(data_surf{ii}.depth_17(end))
    set(f,'CurrentAxes',ha(1))
end
[~, ind_2017] =  min(abs(data_surf{1}.time - datenum(2017,7,1)));
axis tight

% h(5) = errorbar([1 1].*data_surf{1}.time(ind_2017), [6 33.5], [0.1, 2], ...
%     'ok','MarkerFaceColor','k');
h(4) = plot([1 1].*data_surf{1}.time(ind_2017), [32.5 32.5],...
    'ok','MarkerFaceColor','k');
set_monthly_tick(data_surf{4}.time)
% xlim(data_surf{2}.time([ind_2017 end]))
ylim([0 80])
set(gca,'YDir','reverse','XTickLabelRotation',45)
xlabel('Year')
grid on
ylabel('Depth (m)')
d1 = plot(NaN,NaN,'Color',RGB('gray'),'LineWidth',2);
d2 = plot(NaN,NaN,'--','Color',RGB('gray'),'LineWidth',2);
legendflex(h,{model_list2{2:4}},...
    'ref',f,'anchor',{'n' 'n'},...
    'ncol',3,...
    'Location','northoutside','Interpreter','none')
    h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'A');
h_text.FontWeight='bold';
h_text.FontSize=16;

title('Top of the debris layer','FontSize',14)
set(f,'CurrentAxes',ha(2))
title({'Meltwater percolation depth'},'FontSize',14)

hold on
[~, ind_2085] =  min(abs(data_surf{2}.time - datenum(2085,1,1)));
set_monthly_tick(data_surf{4}.time(ind_2085:end))
xlim(data_surf{2}.time([ind_2085 end]))
% ylim([0 80])
set(gca,'YDir','reverse','XTickLabelRotation',45,'YAxisLocation','right')
xlabel('Year')
grid on
ylabel('Depth (m)')

% Create rectangle
h_ann = annotation(f,'rectangle',...
    [ha(1).Position(1)  + ha(1).Position(3)*0.8...
    ha(1).Position(2)+ha(1).Position(4)*0.95 0.07 0.03]);
annotation(f,'line',[h_ann.Position(1)  + h_ann.Position(3)...
    ha(2).Position(1)],...
    [h_ann.Position(2) ha(2).Position(2)]);
annotation(f,'line',[h_ann.Position(1)  + h_ann.Position(3)...
    ha(2).Position(1)],...
    [h_ann.Position(2)+ h_ann.Position(4) ...
    ha(2).Position(2)+ha(2).Position(4)]);
    h_now = gca;
h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
    h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.08, 'B');
h_text.FontWeight='bold';
h_text.FontSize=16;
print(f, sprintf('%s/_10-tracking_debris_horizon',OutputFolder), '-djpeg') 

%% Plotting infiltration for 1966-2100
  f = figure('Visible',vis,'OuterPosition',[0 0 18 14]);
h = [];
ha= tight_subplot(2,1,0.08,[0.3 0.1] ,0.12 );                      
% ha(2).Position(4) = 0.4;
set(f,'CurrentAxes',ha(1))
hold on
for ii =1:4
      
    depth_1 = data_subsurf{ii}.depth_act;
    depth_1(data_subsurf{ii}.slwc<1e-6) = 0;
    perc_depth = max(depth_1,[],1);
    
    plot(data_surf{ii}.time,...
        perc_depth,...
        'Color',col(ii,:),'LineWidth',1.5);
    disp(model_list{ii})
    disp(max(perc_depth))
end
ylabel('perc_depth')
set_monthly_tick(data_surf{ii}.time)
axis tight
% xlim(datenum([2000 2019],1,1))
legendflex(model_list2,'nrow',1,'ref',gcf,'anchor',{'n','n'})

set(f,'CurrentAxes',ha(2))
hold on
for ii =1:4   
    plot(data_surf{ii}.time,...
        data_surf{ii}.melt_mweq*1000,...
        'Color',col(ii,:),'LineWidth',1.5);
    disp(model_list{ii})
    disp(max(perc_depth))
end
ylabel('melt amount')
set_monthly_tick(data_surf{ii}.time)
axis tight
% xlim(datenum([2000 2019],1,1))
%% Plotting 1966 horizon over 1966-2017
figure
hold on
for ii =2:4

       plot(data_surf{ii}.time(1:ind_2017),data_surf{ii}.depth_66(1:ind_2017),...
        'Color',col(ii,:),'LineWidth',2);
       plot(data_surf{ii}.time(1:ind_2017),data_surf{ii}.depth_12(1:ind_2017),...
        'Color',col(ii,:),'LineWidth',2);
end
axis tight
set_monthly_tick(data_surf{ii}.time(1:ind_2017))
%% Comparing SMB with MAR  + Mean melt start end + snowpit
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["year", "smb", "runoff"];
opts.VariableTypes = ["double", "double", "double"];
opts = setvaropts(opts, [1, 2, 3], "DecimalSeparator", ".");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
mar_rcp85 = readtable("Input\Extra\Camp Century\MAR_ACCESS1.3_CC_smb_runoff.txt", opts);
clear opts

disp('-----')
fprintf('Mean SMB (mm w.e.): 2006-2100   2090-2100\n')
fprintf('MAR_ACCESS1.3_RCP8.5 %0.2f  %0.2f\n',...
    mean(mar_rcp85.smb),mean(mar_rcp85.smb(end-10:end)))
for ii = 2:4
    DV = datevec(data_surf_yr{ii}.time);
    fprintf('%s  %0.2f  %0.2f\n',model_list{ii},...
        mean(data_surf_yr{ii}.SMB_mweq(and(DV(:,1)>=2006,...
        DV(:,1)<=2100))*1000),...
        mean(data_surf_yr{ii}.SMB_mweq(and(DV(:,1)>=2090,...
        DV(:,1)<=2100))*1000))
end
disp('-----')

f = figure('vis',vis);
hold on
for ii = 1:4
    DV = datevec(data_surf_yr{ii}.time);
    stairs(DV(:,1),data_surf_yr{ii}.SMB_mweq*1000)
end
stairs(mar_rcp85.year,mar_rcp85.smb,'LineWidth',2)
axis tight; box on; grid on;
set(gca,'XMinorTick','on','YMinorTick','on')
xlabel('Years')
ylabel('Annual surface mass balance (mm w.e.)')
legend({model_list2{:}, 'MAR_ACCESS1.3_RCP8.5'},'Interpreter','none','Location','eastoutside')
print(f, sprintf('%s/SMB_comp_MAR',OutputFolder), '-djpeg')

% Mean melt start end
disp(DV([2, 11]))
disp(DV([end-10 end-1]))
col = lines(4);

disp('model mean_1966-75_mm   std_1966-75_mm     mean_2008-2017_mm   std_2008-2017_mm  mean_2091-2100_mm   std_2091-2100_mm   %change');
for ii = 1:4
    fprintf('%s\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n',model_list{ii},...
        mean(data_surf_yr{ii}.melt_mweq(2:11))*1000000,...
        std(data_surf_yr{ii}.melt_mweq(2:11))*1000000,...
        mean(data_surf_yr{ii}.melt_mweq(43:52))*1000000,...
        std(data_surf_yr{ii}.melt_mweq(43:52))*1000000,...
        mean(data_surf_yr{ii}.melt_mweq((end-10):end-1))*1000000,...
        std(data_surf_yr{ii}.melt_mweq((end-10):end-1))*1000000,...
        (mean(data_surf_yr{ii}.melt_mweq((end-10):end-1)) - mean(data_surf_yr{ii}.melt_mweq(2:11)))/mean(data_surf_yr{ii}.melt_mweq(2:11))*100);
end

% f= figure('OuterPosition',[0 0 13 13]);
% ind_col=[1 1 2 2 3 3 4 4];
% boxplot([tmp],'Color','k',...
%     'labels',{'1966-1976','2007-2017','1966-1976',...
%     '2090-2100','1966-1976','2090-2100','1966-1976',...
%     '2090-2100'})
% set(gca,'XTickLabelRotation',45)
% ylabel('Annual melt (mm w.e.)')
% ha = findall(gca,'Tag','Box');
% for i=1:length(ha)
%     h_patch(i) = patch(ha(i).XData,ha(i).YData,col(ind_col(8-i+1),:));
%     uistack(h_patch(i),'bottom')
% end
% legend(h_patch(8:-2:1), model_list,'Interpreter','none','location','northwest','box','off');
% print(f, sprintf('%s/melt_ini_final',OutputFolder), '-djpeg')

%% Comparison with snow pit data
[pit_data, ~] = ImportSnowpitData(c);

ind = or(strcmp(pit_data.Station,'GITS'),strcmp(pit_data.Station,'CEN'));
pit_data = pit_data(ind,:);
pit_data.date_end = datenum(pit_data.Date);
    temp = datevec(pit_data.date_end);
pit_data.date_start = datenum(temp(:,1)-1,07,15);
pit_data.subl_mod = NaN(size(pit_data,1),4);
pit_data.SMB_mod = NaN(size(pit_data,1),4);

for i=1:size(pit_data,1)
    for ii = 1:4
        % closest time step to the snowpit survey date
        [~, pit_data.ind_end(i)] = min(abs(pit_data.date_end(i) - data_surf{ii}.time));

        % We assume winter accumulation starts on 1st September
        [~, pit_data.ind_start(i)] = min(abs(data_surf{ii}.time-pit_data.date_start(i)));

        % summing the sublimation that occured since the 1st Sept. until
        % the snowpit survey
        pit_data.subl_mod(i,ii) = sum(data_surf{ii}.sublimation_mweq( ...
            pit_data.ind_start(i):pit_data.ind_end(i))); % in m weq

        pit_data.SMB_mod(i,ii) = sum(data_surf{ii}.snowfall(...
            pit_data.ind_start(i):pit_data.ind_end(i))) ...
            + pit_data.subl_mod(i,ii);
    end
end

disp('  ')
disp('Year; Snowpit SWE; RACMO; CanESM_rcp26; CanESM_rcp45; CanESM_rcp85')
tmp = datetime(pit_data.Date);
for i=1:length(tmp.Year)
fprintf(['%0.0f; %0.0f;' num2str(pit_data.SMB_mod(i,:)*1000,'%0.0f; ') '\n'],...
    tmp.Year(i), pit_data.SWE_pit(i))
end
pit_data.SWE_pit(1) = NaN;
fprintf(['ME; ; ' num2str(nanmean(...
    repmat(-pit_data.SWE_pit,1,size(pit_data.SMB_mod,2))...
    + pit_data.SMB_mod*1000,1),'%0.0f; ') '\n' ])
fprintf(['RMSE; ; ' num2str(sqrt(nanmean(...
    (repmat(pit_data.SWE_pit,1,size(pit_data.SMB_mod,2))...
    - pit_data.SMB_mod*1000).^2,1)),'%0.0f; ') '\n' ])

%% Outputting observed subsurface temperature
% disp('Printing corrected observed firn temperature')
% mkdir('.\Output\Corrected observed subsurface temperature')
% clc
% for ii = 1:length(path_list)
%     filename = ['.\Output\Corrected observed subsurface temperature\' ,...
%         station '_T_firn_obs.nc'];
%     disp(station)
%     tic
%     WriteNC_2D(filename, time_mod, data_obs.depth_s,...
%         data_obs.T_ice_obs, 'Depth','m',...
%         'Thermistor depth below the surface', ...
%         'Firn temperature', 'degC', 'Firn temperature'); 
%     toc
% end

%

%% Comparison with GRL2016

opts = spreadsheetImportOptions("NumVariables", 2);
opts.Sheet = "Nor8.5";
opts.DataRange = "A4:B153";
opts.VariableNames = ["year", "smb_obs"];
opts.SelectedVariableNames = ["year", "smb_obs"];
opts.VariableTypes = ["double", "double"];
accum_obs = readtable(".\Input\Extra\Camp Century\MAR_timeseries_calibration_GRL2016.xlsx", opts, "UseExcel", false);
clear opts

opts = spreadsheetImportOptions("NumVariables", 14);
opts.Sheet = "Nor8.5";
opts.DataRange = "A3:N153";
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "year", "SF", "ME", "SMB", "RU", "Var10", "Var11", "SF_adj", "Var13", "SMB_adj"];
opts.SelectedVariableNames = ["year", "SF", "ME", "SMB", "RU", "SF_adj", "SMB_adj"];
opts.VariableTypes = ["char", "char", "char", "char", "double", "double", "double", "double", "double", "char", "char", "double", "char", "double"];
opts = setvaropts(opts, [1, 2, 3, 4, 10, 11, 13], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4, 10, 11, 13], "EmptyFieldRule", "auto");
MAR_NorESM1_RCP85 = readtable(".\Input\Extra\Camp Century\MAR_timeseries_calibration_GRL2016.xlsx", opts, "UseExcel", false);
clear opts

opts = spreadsheetImportOptions("NumVariables", 14);
opts.Sheet = "Can8.5";
opts.DataRange = "A3:N153";
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "year", "SF", "ME", "SMB", "RU", "Var10", "Var11", "SF_adj", "Var13", "SMB_adj"];
opts.SelectedVariableNames = ["year", "SF", "ME", "SMB", "RU", "SF_adj", "SMB_adj"];
opts.VariableTypes = ["char", "char", "char", "char", "double", "double", "double", "double", "double", "char", "char", "double", "char", "double"];
opts = setvaropts(opts, [1, 2, 3, 4, 10, 11, 13], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4, 10, 11, 13], "EmptyFieldRule", "auto");
MAR_CanESM2_RCP85 = readtable(".\Input\Extra\Camp Century\MAR_timeseries_calibration_GRL2016.xlsx", opts, "UseExcel", false);
clear opts

ii = 4;
data_yr = AvgTable(data_surf{ii},'yearly',@sum);
DV = datevec(data_yr.time);
data_yr.year = DV(:,1);
data_yr.SMB = data_yr.snowfall + data_yr.sublimation_mweq;

col =lines(7);%brewermap(7,'lines');
col = col([6 5 4],:);
f = figure('OuterPosition',[0 0 23 11]);
ha = tight_subplot(2,2,0.03,[0.13],[0.1]);
set(f,'CurrentAxes',ha(1))
hold on
h = stairs0(MAR_CanESM2_RCP85.year,MAR_CanESM2_RCP85.SF_adj,'Linewidth',1.5,'Color',col(1,:));
stairs0(MAR_NorESM1_RCP85.year,MAR_NorESM1_RCP85.SF_adj,'Linewidth',1.5,'Color',col(2,:))
stairs0(data_yr.year,data_yr.snowfall*1000,'Linewidth',1.5,'Color',col(3,:))
ylabel({'Snowfall', '(mm w.e.)'})
uistack(h,'top')
box on; axis tight; grid on; set(gca,'XMinorTick','on','XLim',[1966 2100],'FontSize',11)

set(f,'CurrentAxes',ha(2))
hold on
stairs0(MAR_CanESM2_RCP85.year,MAR_CanESM2_RCP85.ME,'Linewidth',1.5,'Color',col(1,:))
stairs0(MAR_NorESM1_RCP85.year,MAR_NorESM1_RCP85.ME,'Linewidth',1.5,'Color',col(2,:))
stairs0(data_yr.year,data_yr.melt_mweq*1000,'Linewidth',1.5,'Color',col(3,:))
stairs0(data_yr.year,data_yr.melt_mweq*1000,'Linewidth',1.5,'Color',col(3,:))
ylabel({'Melt', '(mm w.e.)'})
box on; axis tight; grid on; set(gca,'XMinorTick','on','XLim',[1966 2100],'FontSize',11,'YAxisLocation','right')

set(f,'CurrentAxes',ha(3))
hold on
stairs0(MAR_CanESM2_RCP85.year,MAR_CanESM2_RCP85.RU,'Linewidth',1.5,'Color',col(1,:))
stairs0(MAR_NorESM1_RCP85.year,MAR_NorESM1_RCP85.RU,'Linewidth',1.5,'Color',col(2,:))
stairs0(data_yr.year,data_yr.runoff*1000,'Linewidth',1.5,'Color',col(3,:))
xlabel('Year'); ylabel({'Runoff', '(mm w.e.)'})
box on; axis tight; grid on; set(gca,'XMinorTick','on','XLim',[1966 2100],'FontSize',11)

set(f,'CurrentAxes',ha(4))
hold on
h(1) = stairs0(MAR_CanESM2_RCP85.year,MAR_CanESM2_RCP85.SMB_adj,'Linewidth',1.5,'Color',col(1,:));
h(2) = stairs0(MAR_NorESM1_RCP85.year,MAR_NorESM1_RCP85.SMB_adj,'Linewidth',1.5,'Color',col(2,:));
h(3) = stairs0(data_yr.year,data_yr.SMB_mweq*1000,'Linewidth',1.5,'Color',col(3,:));
uistack(h(1),'top');
xlabel('Year'); ylabel({'Surface mass','balance', '(mm w.e.)'})
box on; axis tight; grid on; set(gca,'XMinorTick','on','XLim',[1966 2100],'FontSize',11,'YAxisLocation','right')

ha(1).XTickLabel = '';
ha(2).XTickLabel = '';

legendflex(h,{'MAR-CanESM2 RCP8.5', 'MAR-NorESM1 RCP8.5', ...
    'GEUS-CanESM2 RCP8.5'},...
    'ref',gcf,...
    'anchor',{'n' 'n'},...
    'box','off',...
    'ncol',3,...
    'interpreter','none');


for i = 1:4
    set(f,'CurrentAxes',ha(i))
    h_now = gca;
    h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
        h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.15, ABC(i));
    h_text.FontWeight='bold';
    h_text.FontSize=10;
end

print(f, sprintf('%s/_6-Comparison_GRL16',OutputFolder), '-djpeg') 


% 
% %     disp(mean(SMB_station(2:42)))
% 
% h(ii) = stairs0(years+month/12,SMB_station);
% h(ii).Color =    col_mod(ii,:);
% h(ii).LineWidth =    2;
% 
% accum_all(accum_all.year<years(1),:) = [];
% fprintf('%i) %s, ME (mm w.e.): %0.2f\n',ii,model_list{ii},...
%     1000*nanmean(SMB_station(1:length(accum_all.median))-accum_all.median));
% fprintf('relative bias (perc): %0.2f\n',...
%     100*nanmean(SMB_station(1:length(accum_all.median))-accum_all.median)/...
%     nanmean(SMB_station(1:length(accum_all.median))));
% fprintf('mean sublimation (mm w.e.): %0.2f\n',...
%     1000*nanmean(data_yr.sublimation_mweq(1:length(accum_all.median))));
% 
%  set(gcf,'CurrentAxes',ha(2))
%  hold on
%     DV = datevec(data_yr.time);
%     years = DV(:,1);
%     SMB_station = data_yr.snowfall + data_yr.sublimation_mweq;
%     h(ii) = stairs0(years(1:end-1),SMB_station(1:end-1));
%     h(ii).Color =    col_mod(ii,:);
%     h(ii).LineWidth =    2;
% end

