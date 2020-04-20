% Script developped by B. Vandecrux (bav@geus.dk) for:
% Vandecrux B, Fausto RS, Langen PL, van As D, MacFerrin M, Colgan WT, 
% Ingeman-Nielsen T, Steffen K, Jensen NS, Møller MT and Box JE (2018) 
% Drivers of Firn Density on the Greenland Ice Sheet Revealed by Weather 
% Station Observations and Modeling. Journal of Geophysical Research: 
% Earth Surface 123(10), 2563–2576. doi:10.1029/2017JF004597.

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
addpath(genpath('lib'))
addpath(genpath('Input'),genpath('Output'))

OutputFolder = './Output/Standard runs';
 folderlist = { './Output/Standard runs/CP1_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10',...
         './Output/Standard runs/DYE-2_long_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10',...
         './Output/Standard runs/NASA-SE_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10',...
        './Output/Standard runs/Summit_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10'};
    vis = 'on';
    station = {'Crawford Point','Dye-2','NASA-SE','Summit'};
    station2 = {'CrawfordPoint','Dye2','NASASE','Summit'};

%% Printing SEB a files if needed
origin_table_1 = table;
origin_table_2 = table;

for i = 1:length(folderlist)
    file_seb_year = sprintf('%s/SEB_year.txt',folderlist{i});
    file_seb_JJA = sprintf('%s/SEB_JJA.txt',folderlist{i});
    
%     if or(~exist(file_seb_year,'file'),~exist(file_seb_JJA,'file'))
        % extract run parameters
        load(strcat(folderlist{i},'/run_param.mat'))
        c.OutputFolder = folderlist{i};
        % extract surface variables
        namefile = sprintf('%s/surf-bin-%i.nc',folderlist{i},1);
        finfo = ncinfo(namefile);
        names={finfo.Variables.Name};
        for ii= 1:size(finfo.Variables,2)
            eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{ii}), namefile,char(names{ii})));
        end
        time_mod = datenum(Year,1,Day,Hour,0,0);
        if ~isfield(c,'verbose')
            c.verbose = 1;
        end
        % extracting observed

        [time_yr, year, day, hour, pres,...
    ~, T, ~, H_instr_temp, o_T1,o_T2, ...
    ~, RH, ~, H_instr_hum, o_RH1, o_RH2, ...
    ~, WS, ~, H_instr_wind, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_AWS, c] = ...
    ExtractAWSData(c);

       SEB_hour = table(time_mod,SHF,LHF,SRin - SRout,LRin - LRout_mdl,...
        rainHF,GF, meltflux,meltflux.*c.dt_obs./c.dev./c.L_fus./c.rho_water,...
        'VariableNames',{'time','SHF_Wm2','LHF_Wm2','SRnet_Wm2','LRnet_Wm2','RainHF_Wm2','GF_Wm2','MeltEnergy_Wm2',...
        'Melt_mweq'});

        SEB_JJA = AvgTableJJA(SEB_hour,'sum');
        SEB_year = AvgTable(SEB_hour,'yearly','sum');
        SEB_year.time = datestr(SEB_year.time);
    
        writetable(SEB_hour, sprintf('%s/SEB_hour.txt', folderlist{i}),'Delimiter',';') 
        writetable(SEB_year, sprintf('%s/SEB_year.txt', folderlist{i}),'Delimiter',';') 
        writetable(SEB_JJA, sprintf('%s/SEB_JJA.txt', folderlist{i}),'Delimiter',';') 
            % Printing tables with origin
            M = [sum(data_AWS.ShortwaveRadiationDownWm2_Origin==0);...
                sum(data_AWS.ShortwaveRadiationUpWm2_Origin==0);...
                sum(data_AWS.AirTemperatureC_Origin==0);...
                sum(data_AWS.RelativeHumidity_Origin==0);...
                sum(data_AWS.AirPressurehPa_Origin==0);...
                sum(data_AWS.WindSpeedms_Origin==0)]./length(data_AWS.WindSpeedms_Origin);
            origin_table_1.(station2{i}) = round(M.*100);
            %     end
end

        writetable(origin_table_1, './Output/origin_table_1.csv','Delimiter',';') 

%% Melt and SEB analysis

clearvars SEB_JJA
for i =1:length(folderlist)
        
    filename = sprintf('%s/SEB_JJA.txt',...
        folderlist{i});
    delimiter = ';';
    startRow = 2;
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    SEB_JJA{i} = table(dataArray{1:end-1}, 'VariableNames', ...
        {'time','SHF','LHF','SRnet','LRnet','RainHF','GF','MeltEnergy','Melt_mweq','Year'});
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
end


years_all = 1998:2015;

FontSize = 16.5;
f = figure('Visible',vis,'units','normalized','outerposition',[0 0 2/3 1]);
ha = tight_subplot(4,2,[0.045 0.01],[0.16 0.12],0.2);

count_plot = 1:2:8;
count_plot_2 = 2:2:8;
for i =1:length(folderlist)

    years_SEB = datevec(SEB_JJA{i}.time);
    years_SEB = years_SEB(:,1);
    SEB_JJA{i}(find(years_SEB>years_all(end)),:) = [];
    SEB_JJA{i}(years_SEB<years_all(1),:) = [];

    set(f,'CurrentAxes',ha(count_plot_2(i)))

    hold on
    x = [];
    temp = datevec(SEB_JJA{i}.time);
    x = temp(:,1);

    x=x';
    disp(station{i})

    fprintf('\tmelt \t%0.2f\t%0.2f\n', nanmean(SEB_JJA{i}.Melt_mweq*1000),...
        nanvar(SEB_JJA{i}.Melt_mweq*1000));
    
    y = NaN(size(years_all));
    ind = ismember(years_all,x);
    y(ind) = SEB_JJA{i}.Melt_mweq*1000;
    x = years_all';
    
%     lm = fitlm(x,y);

    scatter(x,y,40,[0 0 0],'o','fill','LineWidth',2);
    hold on
    plot(x([1 find(~isnan(y),1,'last')]),...
        nanmean(y)*[1 1 ],':k');

box on
set(gca,'XTickLabelRotation',45,'TickLength',[0.02 0.07],'FontSize',FontSize)
axis tight
xlim([x(1)-0.5 x(end)+0.5])
ylim([0 950])
set(gca,'XMinorTick','on','YMinorTick','on')

set (gca,'XTick',x,'XMinorTick','off')
labels=get(gca,'XTickLabel');
labels(1:2:end)={' '};
set(gca,'XTickLabel',labels);
   if ismember(i,1:3)
       set(gca,'XTickLabel','')
   else
       xlabel('Year')
   end
    if i==3
   set(gca,'YAxisLocation','right')

        h_label = ylabel('Melt (mm w.eq.)');
        h_label.Units = 'Normalized';
        h_label.Position = h_label.Position + [0 0.5 0];
   end
   set(gca,'YAxisLocation','right')
end

col = linspecer(6);

hh=[];
symbol = {'' 'o' '^' 'v' 's' '' 'd' };
for i =1:length(folderlist)
    average = [];
    set(f,'CurrentAxes',ha(count_plot(i)))

    hold on
    year_SEB=datetime(datestr(SEB_JJA{i}.time)).Year;
    
    SEB_JJA{i}(year_SEB>max(years_all),:)=[];
    year_SEB=datetime(datestr(SEB_JJA{i}.time)).Year;

    varname=SEB_JJA{i}.Properties.VariableNames;
%     disp(station{i})
    for j = 2:length(varname)-3
       if j==6
        continue
       end
       
        y = NaN(size(years_all));
        ind = ismember(years_all,year_SEB);
        y(ind) = SEB_JJA{i}.(varname{j})./1000000;
%         lm2 = fitlm(SEB_JJA{i}.(varname{j}),SEB_JJA{i}.Melt_mweq);
%     fprintf('\t%s \t%0.2f\t%0.2f\n',varname{j}, nanmean(SEB_JJA{i}.(varname{j})/1000), lm2.Rsquared.Ordinary);

        average = [average, nanmean(SEB_JJA{i}.(varname{j}))];

        lm = fitlm(years_all,y*3600);

         scatter(years_all,y*3600,30,'o','fill',...
            'MarkerEdgeColor',col(j-1,:),'MarkerFaceColor',col(j-1,:),...
            'LineWidth',2);

        fprintf('%s %s %0.2f  (%0.2f) %0.2f %% \n', station{i},varname{j},...
            lm.Coefficients.Estimate(2),...
            max(lm.Coefficients.pValue),...
            (feval(lm,years_all(end))-feval(lm,years_all(1)))/abs(feval(lm,years_all(1)))*100);
% plot(x(ind), ...
%         lm.Coefficients.Estimate(2)*x(ind)+lm.Coefficients.Estimate(1),...
%     'r','LineWidth',2)
    end
    scatter(years_all,NaN*y*3600,50,'o','fill',...
        'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'LineWidth',2);

    hhh = plot([years_all(1) years_all(end)+1],[0 0],'--k');
    if i~=4
    uistack(hhh,'bottom')
    end

    average_pos = average;
    average_pos(average_pos<0) =0;
    
    average_neg = -average;
    average_neg(average_neg<0) =0;
    
set(gca,'XTickLabelRotation',45)

box on
axis tight
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.07],'FontSize',FontSize)
h_title = title(sprintf('%s) %s',...
    char(i+96), station{i}));
h_title.Units = 'Normalized';
h_title.Position = h_title.Position + [0.5 0 0];
h_title.FontSize = FontSize;

set (gca,'XTick',years_all','XTickLabel',num2str(years_all'),'XMinorTick','off')
labels=get(gca,'XTickLabel');
labels(1:2:end,:)=' ';
set(gca,'XTickLabel',labels);
  xlim([years_all(1)-0.5 years_all(end)+0.5])
%   ylim([-30 30])
  ylim([-320 680])

   if ismember(i,1:3)
       set(gca,'XTickLabel','')
   else
       xlabel('Year')
   end
   if ismember(i,[2 4])
%        set(gca,'YTickLabel','')
   elseif i==3
        h_label = ylabel(sprintf('Cumulated JJA energy input (MJ m^{-2})'),'Interpreter','tex');
        h_label.Units = 'Normalized';
        h_label.Position = h_label.Position + [0 0.5 0];
   end
end
legendflex( {'Sensible heat flux','Latent heat flux','Net shortwave radiation',...
    'Net longwave radiation','Conductive heat flux','Melt'}, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 0], ...
                       'nrow',2, ...
                       'fontsize',FontSize,...
                       'box','off');
    uistack(hhh,'bottom')

              
print(f, sprintf('%s/2_SEB_melt',OutputFolder), '-dtiff')
 orient('landscape')
print(f, sprintf('%s/2_SEB_melt',OutputFolder), '-dpdf','-r0')

%% Loading density data
   folders_main = folderlist;
   
    depth_act = cell(3,4);
    rho_all = cell(3,4);
    rho_20 = cell(3,4);
    time_mod = cell(1,4);
    for i = 1:length(folderlist)
        disp(station{i})
        % loading the data from the main run
        [depth_act{1,i}, ~, rho_all{1,i}, time_mod{i}, rho_20{1,i}] = LoadRhoModel(folders_main{i});
        disp('Standard run loaded')   
    end
    
%% plotting density evolution
       f = figure;
    ha = tight_subplot(4,1,0.02,0.07,0.07);
for i =1:length(folderlist)

        set(f,'CurrentAxes',ha(i))
        hold on 

        plot(time_mod{i}, rho_20{1,i},'k','LineWidth',3)

        
        if i ==1
            legendflex( {'Standard run'}, 'ref', gcf, ...
                                   'anchor', {'n','n'}, ...
                                   'buffer',[0 0], ...
                                   'nrow',1, ...
                                   'fontsize',15,...
                                   'box','off');
        end
        time_ext = datenum(1998,7,1):1/24:datenum(2015,6,1);
        set_monthly_tick(time_ext);
        set(gca,'XTickLabelRotation',0)
        if i == 4
            xlabel('Year')
        else
            set(gca,'XTickLabel','')
        end
        if i==3
            h_lab = ylabel('Average density of the upper 20 m of firn (kg/m3)','Interpreter','tex');
            h_lab.Position(2) = h_lab.Position(2)+30;
        end
        axis tight
        ylim([440 650])
          xlim([datenum(1998,7,1) datenum(2015,6,1)]) 

        h_title = title(sprintf('%s) %s',char(i+96),station{i}));
        h_title.HorizontalAlignment = 'left';
        h_title.Position(2) = 600;
        h_title.Position(1) = datenum(2005,1,1);
        set(gca,'layer','top')
    end
    
    print(f,sprintf('%s/Density_evolution',OutputFolder),'-dtiff');

%% Plotting all cores 
    load 'Core_data'

    core_list = table();
    core_list.Site = station';
    core_list.i_core = {FindCore(Core,'Name','G1'):FindCore(Core,'Name','G9');... % CP1
        [FindCore(Core,'Name','core_5_2013') FindCore(Core,'Name','core_6_2013') ...
        FindCore(Core,'Name','core_11_2015')];... %Dye-2
        FindCore(Core,'Name','core_3_2015'):FindCore(Core,'Name','core_6_2015');... % NASA-SE
        [FindCore(Core,'Name','core_22_2015'):FindCore(Core,'Name','core_25_2015') ...
        FindCore(Core,'Name','Albert_2007'):FindCore(Core,'Name','Albert_2000') ];}; %Summit
    
    for i = 1:length(folderlist)
        disp(core_list.Site{i})
        for ii = 1:length (core_list.i_core{i})
            fprintf('    %s %s\n', ...
                Core{core_list.i_core{i}(ii)}.Info.Name, ...
                Core{core_list.i_core{i}(ii)}.Info.NearestCodeLocation)
        end
    end

   for i =1:length(folderlist)    
        % ordering cores in chronological order
        dates = zeros(size(core_list.i_core{i}));
        for ii = 1:length(core_list.i_core{i})
            dates(ii) = datenum(Core{core_list.i_core{i}(ii)}.Info.DateCored);
        end
        [~, i_ordered] = sort(dates);
        core_list.i_core{i} = core_list.i_core{i}(i_ordered);
    
        % removing the cores that were drilled more than 3 years before the
        % beginning of the model run
        i_remove = [];
        for ii = 1:length(core_list.i_core{i})
            temp = datetime(datestr(time_mod{i}(1)));
            if Core{core_list.i_core{i}(ii)}.Info.DateCored.Year < temp.Year - 3
                i_remove = [i_remove, ii];
            end
        end
        core_list.i_core{i}(i_remove) = [];
        num_plot = max(length(core_list.i_core{i}),7);
        f = figure('Visible',vis,'Units','Normalized','outerposition',[0 0 0.9 0.8]);
        [ha, ~] = tight_subplot(1, num_plot, 0.025, [0.15 0.2], [0.09 0.01]);
        count = 0;
                ylim_core = 20;
    
        for ii = core_list.i_core{i}
            count = count+1;
            set(f,'CurrentAxes',ha(count))

            depth = Core{ii}.Data.Depth/100;
            density = Core{ii}.Data.Density;
            h = [];
            hold on
            density2 = density;
            depth2=depth;
            ind_last = find(~isnan(density2),1,'last');
            density2=density2(1:ind_last);
            depth2=depth2(1:ind_last);
            ind_first = find(~isnan(density2),1,'first');
            density2= density2(ind_first:end);
            depth2= depth2(ind_first:end);
            if isnan(density(1))
                density2(1) = 315;
            end
            density2(isnan(density2))=interp1(find(~isnan(density2)),...
                density2(~isnan(density2)),...
                find(isnan(density2)),'linear');
            if isrow(density2)
                density2=density2';
            end
            if isrow(depth2)
                depth2=depth2';
            end
            
    x2 = [density2-23.1*2; flipud(density2+23.1*2)];
    inBetween = [depth2; flipud(depth2)];

    
    fill(x2, inBetween, 0.8*[1 1 0] + [0 0 1],'LineStyle','none');
            h(1) = plot(density,depth,'b','LineWidth',2);
            
            time_core = datenum(Core{ii}.Info.DateCored);
            [~, ind_time] = min(abs( time_mod{i} - time_core));
        
            depth_mod = depth_act{1,i}(:,ind_time);
            density_mod = rho_all{1,i}(:,ind_time);
            h(2) = stairs([density_mod; density_mod(end)], [0; depth_mod],...
                'k', 'LineWidth',2);
        
            set(gca,'Ydir','reverse','XMinorTick','on','YMinorTick','on','XTickLabelRotation',45)
            xlim([100 900])
            set(gca,'XTick',[100 500 900])
            ylim([0 ylim_core ])
            
            h_title = title(sprintf('%s\n%i',Core{ii}.Info.Name,Core{ii}.Info.DateCored.Year));
            h_title.FontSize = 13;
            box on
            if count == 1  
                h_legend = legend(h,'Observed','Modelled');
                set(h_legend,'Parent',gcf)
                h_legend.Units = 'Normalized';
                h_legend.Box ='off';
                h_legend.Position = [0.8*length(core_list.i_core{i})/num_plot 0.95 0.03 0];                
                
                ylabel('Depth (m)')
                AxesH = axes('Parent', gcf, 'Units','Normalized' ,...
                    'Position', [0.02, 0.02, 0.98, 0.98], ...
                    'Visible','off');
                h_text = text(0,1, station{i}, ...
                  'HorizontalAlignment', 'center', ...
                  'VerticalAlignment', 'top');
                h_text.FontSize = 15;
                h_text.FontWeight = 'bold';
                h_text.Units = 'Normalized';
                h_text.Position = [0.5*length(core_list.i_core{i})/num_plot 0.97 0];
                
                h_text = text(0,1, sprintf('Density (kg m^{-3})'), ...
                  'HorizontalAlignment', 'center', ...
                  'VerticalAlignment', 'top',...
                  'Interpreter','tex');
                h_text.FontSize = 15;
                h_text.Units = 'Normalized';
                h_text.Position = [0.5*length(core_list.i_core{i})/num_plot 0.04 0];
                
                           else
                set(gca,'YTickLabel',' ')
            end
            
        end
        for ii = (count+1):num_plot
            set(f,'CurrentAxes',ha(ii))
            set(gca,'Visible','off')
        end

        orient('landscape')
        print(f,...
            sprintf('%s/S8_%s', ...
            OutputFolder,station{i}),'-dtiff');
        print(f,...
            sprintf('%s/S8_%s', ...
            OutputFolder,station{i}),'-dpdf','-r0');
   end
   
%% Scatter plot modelled observed

    shape = {'o', 'd','^','v'};
all_up_meas = [];
all_up_mod = [];
all_down_meas = [];
all_down_mod = [];
h = [];
diff = [];
diff_sqrd = [];

f = figure('Visible',vis);
hold on
    for i= 1:length(folderlist)
        core_list.upper_density_meas{i} = [];
        core_list.lower_density_meas{i} = [];
        core_list.upper_density_mod{i} = [];
        core_list.lower_density_mod{i} = [];
        for ii = core_list.i_core{i}

            % extracting density profile from model
            time_core = datenum(Core{ii}.Info.DateCored);
            [~, ind_time] = min(abs( time_mod{i} - time_core));
            depth_mod = depth_act{1,i}(:,ind_time); % we round at the closest cm
            density_mod = rho_all{1,i}(:,ind_time);
            
            new_depth_mod = 1:floor(depth_mod(end)*100);
            new_density_mod = NaN(size(Core{ii}.Data.Depth));
        
        for jk = length(depth_mod):-1:1
                ind_depth = new_depth_mod <= depth_mod(jk)*100;
                new_density_mod(ind_depth) = density_mod(jk);
        end
        new_depth_mod(new_depth_mod>Core{ii}.Data.Depth(end)) = [];
        new_density_mod(new_depth_mod>Core{ii}.Data.Depth(end)) = [];
           
            % now finding the indexes in the core
            ind_up = and(~isnan(Core{ii}.Data.Density), ...
                Core{ii}.Data.Depth<=500);
            ind_down = and(and(~isnan(Core{ii}.Data.Density),...
                Core{ii}.Data.Depth>500),...
                Core{ii}.Data.Depth <= new_depth_mod(end)*100);
            ind_all = or(ind_up,ind_down);
            diff = [diff mean(new_density_mod(ind_all)) - mean(Core{ii}.Data.Density(ind_all))];
            diff_sqrd = [diff_sqrd (mean(new_density_mod(ind_all)) - mean(Core{ii}.Data.Density(ind_all)))^2];
            
            core_list.upper_density_mod{i} = ...
                [core_list.upper_density_mod{i} mean(new_density_mod(ind_up))];
            core_list.lower_density_mod{i} = ...
                [core_list.lower_density_mod{i} mean(new_density_mod(ind_down))];

            core_list.upper_density_meas{i} = ...
                [core_list.upper_density_meas{i} mean(Core{ii}.Data.Density(ind_up))];
            core_list.lower_density_meas{i} = ...
                [core_list.lower_density_meas{i} mean(Core{ii}.Data.Density(ind_down))];
            
        end


         err = 40*ones(size(core_list.upper_density_meas{i}));
         e=herrorbar(core_list.upper_density_meas{i}, core_list.upper_density_mod{i}, err);
           e(2).LineStyle = 'none';
            e(1).Color = 'k';

         err = 40*ones(size(core_list.lower_density_meas{i}));
         e = herrorbar(core_list.lower_density_meas{i}, core_list.lower_density_mod{i}, err);
           e(2).LineStyle = 'none';
            e(1).Color = 'k';
        h((i-1)*2 + 1) = scatter(core_list.upper_density_meas{i}, core_list.upper_density_mod{i},...
            120,shape{i},'filled','MarkerEdgeColor','b','MarkerFaceColor','b');
        h((i-1)*2 + 2) = scatter(core_list.lower_density_meas{i}, core_list.lower_density_mod{i},...
            120,shape{i},'filled','MarkerEdgeColor','r','MarkerFaceColor','r');
        all_up_meas = [all_up_meas core_list.upper_density_meas{i}];
        all_up_mod = [all_up_mod core_list.upper_density_mod{i}];
        all_down_meas = [all_down_meas core_list.lower_density_meas{i}];
        all_down_mod = [all_down_mod core_list.lower_density_mod{i}];
    end

        plot([280 710], [280 710], 'k')
        axis tight square
        box on
        set(gca,'XMinorTick','on','YMinorTick','on')
        xlabel(sprintf('Measured density (kg m^{-3})'),'Interpreter','tex')
        ylabel(sprintf('Modelled density (kg m^{-3})'),'Interpreter','tex')
    legend(h,'Crawford Point (< 5 m)',...
        'Crawford Point (> 5 m)',...
        'Dye-2 (< 5 m)',...
        'Dye-2 (> 5 m)',...
        'NASA-SE (< 5 m)',...
        'NASA-SE (> 5 m)',...
        'Summit (< 5 m)',...
        'Summit (> 5 m)',...
        'Location','EastOutside')
    title(sprintf('Overall: ME = %0.1f  RMSE = %0.1f\nUpper firn: ME = %0.1f  RMSE = %0.1f\n Deeper firn: ME = %0.1f  RMSE = %0.1f',...
            nanmean(diff), sqrt(nanmean(diff_sqrd)),...
            nanmean(all_up_mod-all_up_meas), sqrt(nanmean((all_up_mod-all_up_meas).^2)),...
            nanmean(all_down_mod-all_down_meas), sqrt(nanmean((all_down_mod - all_down_meas).^2))));

        orient('landscape')
        print(f, sprintf('%s/4_ScatterDensity',OutputFolder), '-dtiff')
        print(f, sprintf('%s/4_ScatterDensity',OutputFolder), '-dpdf','-r0')
        
%% start and end density
disp('Site Start_time End_time Start_rho End_rho Difference  Start_PV End_PV Difference')
for i = 1:length(folderlist)
    ps = @(m,rho) max(0,m.*(1./rho -1/873));

    text = strcat( [station{i}, ('\t'),...
                datestr(time_mod{i}(1)), ('\t'),...
                datestr(time_mod{i}(end)), ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f (%0.2f)', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f (%0.2f)', '\n']);
   fprintf(text,...
       rho_20{1,i}(1),...
        rho_20{1,i}(end), ...
        rho_20{1,i}(end)-rho_20{1,i}(1), ...
(rho_20{1,i}(end)-rho_20{1,i}(1))/rho_20{1,i}(1)*100, ...
        ps(rho_20{1,i}(1)*20,rho_20{1,i}(1)), ...
        ps(rho_20{1,i}(end)*20,rho_20{1,i}(end)), ...
        ps(rho_20{1,i}(end)*20,rho_20{1,i}(end)) - ...
        ps(rho_20{1,i}(1)*20,rho_20{1,i}(1)),...
        (ps(rho_20{1,i}(end)*20,rho_20{1,i}(end)) - ...
        ps(rho_20{1,i}(1)*20,rho_20{1,i}(1)))/ ...
        ps(rho_20{1,i}(1)*20,rho_20{1,i}(1))*100);
end
disp('')
disp('================================')
disp('Site Start_time End_time Start_rho End_rho Difference  Start_PV End_PV Difference')

 for i = 1:2
    ps = @(m,rho) max(0,m.*(1./rho -1/873));

    ind = find(time_mod{i}==datenum(2005,6,1));
    text = strcat( [station{i}, ('\t'),...
                datestr(time_mod{i}(1)), ('\t'),...
                datestr(time_mod{i}(ind)), ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f (%0.2f)', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f (%0.2f)', '\n']);
   fprintf(text,...
       rho_20{1,i}(1),...
        rho_20{1,i}(ind), ...
        rho_20{1,i}(ind)-rho_20{1,i}(1), ...
(rho_20{1,i}(ind)-rho_20{1,i}(1))/rho_20{1,i}(1)*100, ...
        ps(rho_20{1,i}(1)*20,rho_20{1,i}(1)), ...
        ps(rho_20{1,i}(ind)*20,rho_20{1,i}(ind)), ...
        ps(rho_20{1,i}(ind)*20,rho_20{1,i}(ind)) - ...
        ps(rho_20{1,i}(1)*20,rho_20{1,i}(1)),...
        (ps(rho_20{1,i}(ind)*20,rho_20{1,i}(ind)) - ...
        ps(rho_20{1,i}(1)*20,rho_20{1,i}(1)))/ ...
        ps(rho_20{1,i}(1)*20,rho_20{1,i}(1))*100);
    
        text = strcat( [station{i}, ('\t'),...
                datestr(time_mod{i}(ind)), ('\t'),...
                datestr(time_mod{i}(end)), ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f (%0.2f)', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f', ('\t'),...
                '%0.2f (%0.2f)', '\n']);
   fprintf(text,...
       rho_20{1,i}(ind),...
        rho_20{1,i}(end), ...
        rho_20{1,i}(end)-rho_20{1,i}(ind), ...
(rho_20{1,i}(end)-rho_20{1,i}(ind))/rho_20{1,i}(ind)*100, ...
        ps(rho_20{1,i}(ind)*20,rho_20{1,i}(ind)), ...
        ps(rho_20{1,i}(end)*20,rho_20{1,i}(end)), ...
        ps(rho_20{1,i}(end)*20,rho_20{1,i}(end)) - ...
        ps(rho_20{1,i}(ind)*20,rho_20{1,i}(ind)),...
        (ps(rho_20{1,i}(end)*20,rho_20{1,i}(end)) - ...
        ps(rho_20{1,i}(ind)*20,rho_20{1,i}(ind)))/ ...
        ps(rho_20{1,i}(ind)*20,rho_20{1,i}(ind))*100);
 end
 
 %% Density comparison t1 t2
for i = 1:length(folderlist)
    % comparison with the core period
    disp(station{i})
    switch station{i}
        case 'Crawford Point'
            t1 = datenum(1998,6,1);
            t2 = datenum(2007,6,1);
            
            t1 = datenum(1998,6,1);
            t2 = datenum(2007,6,1);
            
            % observed change
            i1 = FindCore(Core,'Name','CORE 6945');           
            i2 = FindCore(Core,'Name','G1'):FindCore(Core,'Name','G9');
            i3 =[];
            i4=[];
        case 'Dye-2'
            t1(1) = datenum(1998,6,1);
            t2(1) = datenum(2013,5,1);
            % observed change
            i1 = FindCore(Core,'Name','DYE2 1998 core B');           
            i2 = FindCore(Core,'Name','core_5_2013'):FindCore(Core,'Name','core_6_2013');
            
            t1(2) = datenum(2013,5,1);
            t2(2) = datenum(2015,5,1);
            % observed change
            i3 = i2;
            i4 = FindCore(Core,'Name','core_11_2015');
        case 'NASA-SE'
            t1 = datenum(1998,6,1);
            t2 = datenum(2015,5,1);
            % observed change
            i1 = FindCore(Core,'Name','CORE 6642 (B)');           
            i2 = FindCore(Core,'Name','core_3_2015'):FindCore(Core,'Name','core_6_2015');
            i3 =[];
            i4=[];
        case 'Summit'
            t1(1) = datenum(2000,7,1);
            t2(1) = datenum(2007,7,1);
            % observed change
            if ~isempty(FindCore(Core,'Name','Dome GRIP'))
                Core{FindCore(Core,'Name','Dome GRIP')}.Data.Density(1:453) = ...
                    Core{FindCore(Core,'Name','Mayewski_1990')}.Data.Density(1:453);
                Core{FindCore(Core,'Name','Dome GRIP')}.Info.Name = 'GRIP & Mayewski';
            end
            i1 = FindCore(Core,'Name','GRIP & Mayewski');           
            i2 = FindCore(Core,'Name','Albert_2007');
            
            t1(2) = datenum(2007,7,1);
            t2(2) = datenum(2015,6,1);
            % observed change
            i3 = i2;           
            i4 = FindCore(Core,'Name','core_22_2015'):FindCore(Core,'Name','core_25_2015');
    end

    % Observed change in density
    depth_range = [];
    for ii=1:2
        if ii == 2
            if isempty(i3)
                continue
            else
                i1 = i3;
                i2 = i4;
            end
        end
        
        ind_last = 0;
        if length(i1)>1
            for kk = i1
                ind_last = max(ind_last, ...
                    find(~isnan(Core{kk}.Data.Density),1,'last'));
            end

            core_1.depth = 1:ind_last;
            core_1.density = NaN * core_1.depth;
            mat = NaN(length(core_1.density),length(i1));
            for j = i1
                ind_max = min(length(core_1.density), ...
                    length(Core{j}.Data.Depth));
                mat(1:ind_max,j-i1(1)+1) = Core{j}.Data.Density(1:ind_max)';
            end
            core_1.density = nanmean(mat,2);
            
        else
            core_1.depth = Core{i1}.Data.Depth;
            core_1.density = Core{i1}.Data.Density;
            if isrow(core_1.depth)
                core_1.depth=core_1.depth';
            end
            if isrow(core_1.density)
                core_1.density=core_1.density';
            end
        end
        core_1.density(core_1.depth>2000) = [];
        core_1.depth(core_1.depth>2000) = [];
        
        core_2.depth = 1:length(core_1.depth)';
        core_2.density = NaN * core_2.depth;
        mat = NaN(length(core_2.density),length(i2));
        for j = i2
            ind_max = min(length(core_2.density), ...
                length(Core{j}.Data.Depth));
            mat(1:ind_max,j-i2(1)+1) = Core{j}.Data.Density(1:ind_max)';
        end
        core_2.density = nanmean(mat,2);

        ind_last = find(~isnan(core_2.density),1,'last');
        core_1.density=core_1.density (1:ind_last);
        core_2.density=core_2.density (1:ind_last);

        %gap-filling
        ind_nan = find(isnan(core_2.density));
        ind_nonan = find(~isnan(core_2.density));
        core_2.density(ind_nan) = ...
            interp1([0; ind_nonan],[315; core_2.density(ind_nonan)],ind_nan,'linear');

        ind_nan = find(isnan(core_1.density));
        ind_nonan = find(~isnan(core_1.density));
        core_1.density(ind_nan) = ...
            interp1([0; ind_nonan],[315; core_1.density(ind_nonan)],ind_nan,'linear');

        disp('Year 1  Year 2  Difference')
        fprintf('%i\t%i\t%0.2f\n',...
            Core{i1(1)}.Info.DateCored.Year,...
            Core{i2(1)}.Info.DateCored.Year,...
            mean(core_2.density) - mean(core_1.density))

        figure
        subplot(2,1,1)
        plot(core_1.density,'LineWidth',2)
        hold on
        for k = i1
            plot(Core{k}.Data.Density)
        end
        xlim([0 2000])
        subplot(2,1,2)
        plot(core_2.density,'LineWidth',2)
        hold on
        for k = i2
            plot(Core{k}.Data.Density)
        end
        xlim([0 2000])
        
        depth_range = [depth_range length(core_2.density)];
    end

    
    for j = 1:length(t1)
        [~, ind_start] = min(abs(time_mod{i}-t1(j)));
        [~, ind_end] = min(abs(time_mod{i}-t2(j)));

            % extracting density profile from model
            % start time
        depth_mod_start = depth_act{1,i}(:,ind_start); % we round at the closest cm
        density_mod_start = rho_all{1,i}(:,ind_start);

        new_depth_mod_start = 1:floor(depth_mod_start(end)*100);
        new_density_mod_start = NaN(size(Core{ii}.Data.Depth));

        for jk = length(depth_mod_start):-1:1
                ind_depth = new_depth_mod_start <= depth_mod_start(jk)*100;
                new_density_mod_start(ind_depth) = density_mod_start(jk);
        end

        rho_avg_1 = mean(new_density_mod_start(new_depth_mod_start<=2000));
            % end time
        depth_mod_end = depth_act{1,i}(:,ind_end); % we round at the closest cm
        density_mod_end = rho_all{1,i}(:,ind_end);

        new_depth_mod_end = 1:floor(depth_mod_end(end)*100);
        new_density_mod_end = NaN(size(Core{ii}.Data.Depth));

        for jk = length(depth_mod_end):-1:1
                ind_depth = new_depth_mod_end <= depth_mod_end(jk)*100;
                new_density_mod_end(ind_depth) = density_mod_end(jk);
        end

        rho_avg_2 = mean(new_density_mod_end(new_depth_mod_end<=2000));


        f = figure('Visible',vis);
        hold on
        plot(new_depth_mod_start'./100, new_density_mod_start,'r','Linewidth',2)
        plot([0 20],[rho_avg_1 rho_avg_1],'--r','Linewidth',2)
        plot(new_depth_mod_end'./100, new_density_mod_end,'b','Linewidth',2)
        plot([0 20],[rho_avg_2 rho_avg_2],'--b','Linewidth',2)

        title(sprintf('%s\n \\Delta\\rho = %0.2f',...
            station{i}, rho_avg_2-rho_avg_1),'Interpreter','tex')
        
        fprintf('Model (20 m): %0.2f\n',rho_avg_2-rho_avg_1)

        rho_avg_lim_1 = mean(new_density_mod_start(new_depth_mod_start<=depth_range(j)));
        rho_avg_lim_2 = mean(new_density_mod_end(new_depth_mod_end<=depth_range(j)));

        fprintf('Model (0 - %0.1f m): %0.2f\n\n',depth_range(j)/100,rho_avg_lim_2-rho_avg_lim_1)

        axis tight square
        xlim([0 20])
        view([90 90])
        xlabel('Depth (m)')
        ylabel('Density (kg m^{-3})')
        box on
        set(gca,'XMinorTick','on','XTick',0:20,...
            'XTickLabel',{0, '', 2, '', 4, '', 6, '', 8, '', 10, '', 12, '', 14 '',16,'', 18,'',20},...
            'fontsize',18)
        legend(datestr(time_mod{i}(ind_start),'dd-mmm-yyyy'),sprintf('mean upper 20 m: %0.2f kg m^{-3}',rho_avg_1),...
            datestr(time_mod{i}(ind_end),'dd-mmm-yyyy'),sprintf('mean upper 20 m: %0.2f kg m^{-3}',rho_avg_2),'Location','EastOutside')
        print(f, sprintf('%s/firn_evolution_%s',OutputFolder,station{i}), '-dtiff')
    end
end
                 
 %% Density plot
  f=figure;%('outerposition',[1 -1  25 25]);
 ha = tight_subplot(4,1,0.02,[0.08 0.02],[0.07 0.07]);
ylimit = 20;
 for ii =1:length(folderlist)
      % extract run parameters
    load(strcat(folderlist{ii},'/run_param.mat'))

    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',folderlist{ii},1);
    Year = ncread(namefile,'Year');
    Day = ncread(namefile,'Day');
    Hour = ncread(namefile,'Hour');

    % extract subsurface variables
    namefile = sprintf('%s/subsurf-bin-%i.nc',folderlist{ii},1);       
    snowc = ncread(namefile,'snowc');
    snic = ncread(namefile,'snic');
    slwc = ncread(namefile,'slwc');
    rho = ncread(namefile,'rho');

    % Update BV2017: (not used anymore) Liquid water does not participate to the actual thickness
    % of a layer. Ice does not participate as long as there is enough pore
    % space in the snow to accomodate it.
    pore_space = snowc .* c.rho_water.*( 1./rho - 1/c.rho_ice);
    excess_ice = max(0, snic * c.rho_water / c.rho_ice - pore_space);
    thickness_act = snowc.*(c.rho_water./rho) + excess_ice;

    % thickness_act = snowc.*(c.rho_water./rho) + snic.*(c.rho_water./c.rho_ice);

    depth_act=zeros(size(thickness_act));
    for j=1:length(thickness_act(1,:))
            depth_act(1:end,j)=cumsum(thickness_act(1:end,j));
    end

    rho_all= (snowc + snic)./...
                (snowc./rho + snic./c.rho_ice);
    lwc = slwc(:,:) ./ thickness_act(:,:);

    time_mod = datenum(Year,1,Day,Hour,0,0);
    TT = ones(c.jpgrnd,1) * time_mod';

    H_surf = depth_act(end,:)'-depth_act(end,1); %+snowbkt*1000/315;

    for i = 1:length(H_surf)-1
        if (H_surf(i+1)-H_surf(i))> c.new_bottom_lay-1
            H_surf(i+1:end) = H_surf(i+1:end) - c.new_bottom_lay*c.rho_water/c.rho_ice;
        end
    end    

    depth_act = vertcat(zeros(size(depth_act(1,:))), depth_act);

    for i = 1:size(depth_act,1)
        depth_act(i,:) = depth_act(i,:) - H_surf' + H_surf(1);
    end


    [~, ind_bot] = min(abs(depth_act - (ylimit+1)));
    ind_lim = max(ind_bot);
    step = 72;
    set(f,'CurrentAxes',ha(ii))
    col = PlotTemp(TT(1:ind_lim, 1:step:end),...:), ...
        depth_act(1:ind_lim, 1:step:end),...:), ...
        rho_all(1:ind_lim, 1:step:end),...:), ...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'yes',...
        'ValueIsoTherm', c.rho_pco,...
        'ShowLegend','no',...
        'cmap','hsv2',...
        ...'Interp','on',...
        'XLabel','',...
        'YLabel','Depth (m)',...
        'CLabel','Firn density (kg m^{-3})',...
        'Range', 300:5:900);
%     col.FontSize = 12;
    if ii == 1
        col.Position(1) = col.Position(1)- 0.03;
        col.Position(4) = 0.9 ;
        col.Position(2) = 0.08;
    else
        col.Position(1) = col.Position(1)+ 3;
    end
    plot(time_mod,-H_surf+H_surf(1),'LineWidth',2)
    ylim([min([-H_surf+H_surf(1)])-1 ylimit]) %max(max(depth_act))
    col.TickLength=col.TickLength*5;
    temp = col.YTickLabel;
    for k = 1:length(temp)
        if (k-1)/10==floor((k-1)/10)
        else
            temp(k,:)=' ';
        end
    end
    col.YTickLabel=temp;

    xlim([datenum(1998,6,1) datenum(2015,6,1)]) 
    time_ext = datenum(1998,6,1):1/24:datenum(2015,6,1);
    set_monthly_tick(time_ext);
    temp = get(gca,'XTickLabel');
    for k = 1:length(temp)
        if k/2==floor(k/2)
            temp(k,:)=' ';
        end
    end
    set(gca,'XTickLabel',temp,'XTickLabelRotation',0);

    if ii~=4
        set(gca,'XTickLabel','')
    end
    clearvars text
    h_text = text(time_ext(24*60), ...
        min([-H_surf+H_surf(1)])-1+abs(min([-H_surf+H_surf(1)])-1)/2,...
        sprintf('%s) %s',char(ii+96),station{ii}));
    h_text.FontSize = 15;
    h_text.FontWeight = 'bold';
%     h_text.Parent = gcf;
%     h_text.Units = 'normalized';
%     h_text.Position = 
 end
 xlabel('Year')
print(f, sprintf('%s/3_density_evolution',OutputFolder), '-dtiff')
% print(f, sprintf('%s/rho_comp2',OutputFolder), '-dpdf')

%% Densification yearly plot all sites
  f=figure;%('outerposition',[1 -1  25 25]);
 ha = tight_subplot(4,1,0.04,[0.08 0.15],[0.1 0.07]);
vis = 'on';
 for ii =1:length(folderlist)
    load(strcat(folderlist{ii},'/run_param.mat'))
    c.OutputFolder = folderlist{ii};
    % extract surface variables
    namefile = sprintf('%s/track_density.nc',folderlist{ii});
    finfo = ncinfo(namefile);
    names={finfo.Variables.Name};
    for i= 1:length(names)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end
    
      % extract run parameters
    load(strcat(folderlist{ii},'/run_param.mat'))
    c.OutputFolder = folderlist{ii};
    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',folderlist{ii},1);
    names = {'Year' 'Day' 'Hour' 'H_surf' 'SRout_mdl' 'LRout_mdl' 'SHF' 'LHF' ...
         'GF' 'rainHF' 'meltflux'   'runoff' 'snowfall' 'rainfall'...
          'sublimation' 'snowbkt'};
    for i= 1:length(names)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end

time_mod = datenum(Year,1,Day,Hour,0,0);
[~, year, day, hour, pres,...
    T1, T2, z_T1, z_T2, o_T1,o_T2, ...
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, ...
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_out, c] = ExtractAWSData(c);


    runoffhour = -[runoff(1); runoff(2:end)-runoff(1:end-1)] ;
    time = time_mod;    
    set(f,'CurrentAxes',ha(ii))

    disp('Plotting densification drivers')
    c.plot_month_dens = 0;

    [h_legend] = DensificationStudy(time_mod,density_avg_20, ...
                                 vis, c);

    xlim([datenum(1998,6,1) datenum(2015,6,1)]) 
    time_ext = datenum(1998,6,1):1/24:datenum(2015,6,1);
    set_monthly_tick(time_ext);
    temp = get(gca,'XTickLabel');
    for k = 1:length(temp)
        if k/2==floor(k/2)
            temp(k,:)=' ';
        end
    end
    set(gca,'XTickLabel',temp,'XTickLabelRotation',0);
    if ii==1
        h_ylab = ylabel('Contribution to density change (kg m^{-3})','Interpreter','tex');

        h_ylab.Units = 'normalized';
        h_ylab.Position(2) =- 1.2;
        h_ylab.Position(1) = -0.06;
    else
        cla(h_legend);
        ylabel('')
    end
    if ii~=4
        set(gca,'XTickLabel','')
        xlabel('')
    end
    ylim([-50 50])
    clearvars text
    h_text = text(time_ext(24*60), ...
        60,...
        sprintf('%s) %s',char(ii+96),station{ii}));
    h_text.FontSize = 15;
    h_text.FontWeight = 'bold';
    h_text.Position(2) = 60;
 end
 
 xlabel('Year')
 set(gcf, 'Units','centimeter','outerposition', [.25 .25 [25 32]-0.5]);
set(gcf,'papersize',[25 32 ]);
 orient('landscape')
print(f, sprintf('%s/5_densification_comp',OutputFolder), '-dtiff')
print(f, sprintf('%s/5_densification_comp',OutputFolder), '-dpdf')

%% Densification hourly plot all sites
  f=figure;%('outerposition',[1 -1  25 25]);
 ha = tight_subplot(4,1,0.04,[0.1 0.08],[0.09 0.07]);
vis = 'on';
 for ii =1:length(folderlist)
    load(strcat(folderlist{ii},'/run_param.mat'))
    c.OutputFolder = folderlist{ii};
    % extract surface variables
    namefile = sprintf('%s/track_density.nc',folderlist{ii});
    finfo = ncinfo(namefile);
    names={finfo.Variables.Name};
    for i= 1:length(names)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end
    
      % extract run parameters
    load(strcat(folderlist{ii},'/run_param.mat'))
    c.OutputFolder = folderlist{ii};
    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',folderlist{ii},1);
    %             finfo = ncinfo(namefile);
    %             names={finfo.Variables.Name};
    names = {'Year' 'Day' 'Hour' 'H_surf' 'SRout_mdl' 'LRout_mdl' 'SHF' 'LHF' ...
         'GF' 'rainHF' 'meltflux'   'runoff' 'snowfall' 'rainfall'...
          'sublimation' 'snowbkt'};
    for i= 1:length(names)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end

time_mod = datenum(Year,1,Day,Hour,0,0);
[time, year, day, hour, pres,...
    T1, T2, z_T1, z_T2, o_T1,o_T2, ...
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, ...
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_out, c] = ExtractAWSData(c);


    runoffhour = -[runoff(1); runoff(2:end)-runoff(1:end-1)] ;
    SMB = snowfall+rainfall+sublimation+runoffhour;
    time = time_mod;
    SMB_hour = table(time, snowfall, rainfall, sublimation, runoffhour,SMB);
    SMB_hour.Properties.VariableNames{5} = 'runoff';
    SMB_daily = AvgTable(SMB_hour,'daily');

    SEB_hour = table(time_mod,SHF,LHF,SRin - SRout,LRin - LRout_mdl,...
        rainHF,GF, meltflux,meltflux.*c.dt_obs./c.dev./c.L_fus./c.rho_water,...
        'VariableNames',{'time','SHF','LHF','SRnet','LRnet','RainHF','GF','MeltEnergy',...
        'Melt_mweq'});

    SEB_daily = AvgTable(SEB_hour,'daily');
            
    set(f,'CurrentAxes',ha(ii))

    disp('Plotting densification drivers')
    c.plot_month_dens = 0;

    rho_avg = density_avg_20(6,:);

    [ax,H1,H2] = plotyy(time_mod,rho_avg,...
        SMB_daily.time,(SMB_daily.snowfall + max(0,SMB_daily.sublimation))*1000);
    H2.LineWidth = 2;
    H2.Color = 'b';
    H1.LineWidth = 4;
    H1.Color = 'k';
    set(ax,{'ycolor'},{'k';'k'}) 
    set(f,'CurrentAxes',ax(1));
    hold on
%     x2 = [time_mod; flipud(time_mod)];
%     inBetween = [rho_avg'+40; flipud(rho_avg')-40];
%     h = fill(x2, inBetween, 0.8*[1 1 1]);
    plot(time_mod,rho_avg,'k','Linewidth',4)
    axis tight
    ylim([460 640])
    set(gca,'YTick', 460:60:640)
    time_ext = datenum(1998,6,1):1/24:datenum(2015,6,1);
    set_monthly_tick(time_ext);
    set(gca,'XTickLabelRotation',0)

        xlim([datenum(1998,6,1) datenum(2015,6,1)]) 
    if ii~=4
        set(gca,'XTickLabel','')
        xlabel('')
    end
        if ii==1
        h_ylab = ylabel('Firn density (kg m^{-3})','Interpreter','tex');

        h_ylab.Units = 'normalized';
        h_ylab.Position(2) =- 1.2;
        h_ylab.Position(1) = - 0.06 ;
    else
        ylabel('')
    end
    
    
    set(f,'CurrentAxes',ax(2));
    hold on
    h3 = plot(SEB_daily.time,SEB_daily.Melt_mweq*1000,'r','LineWidth',2);
    axis tight   
    ylim([0 45])

    xlim([datenum(1998,6,1) datenum(2015,6,1)]) 
    set_monthly_tick(time_ext);
    set(gca,'YTick',0:15:45)
    temp = get(gca,'XTickLabel');
    for k = 1:length(temp)
        if k/2==floor(k/2)
            temp(k,:)=' ';
        end
    end
    set(gca,'XTickLabel',temp,'XTickLabelRotation',0);
    if ii==1
        h_ylab = ylabel('Daily snowfall & melt (mm)');

        h_ylab.Units = 'normalized';
        h_ylab.Position(2) =- 1.2;
legendflex([H1 H2 h3], {'Average 20 m firn density','Snowfall','Melt'}, ...
                        'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 0], ...
                       'nrow',1, ...
                       'fontsize',15,...
                       'box','off');
    else
        ylabel('')
    end
    if ii~=4
        set(gca,'XTickLabel','')
        xlabel('')
    else
                xlabel('Year')

    end
%     ylim([-70 70])
    h_text = text(time_ext(24*60), ...
        39,...
        sprintf('%s) %s',char(ii+96),station{ii}));
    h_text.FontSize = 15;
    h_text.FontWeight = 'bold';
    h_text.Position(2) = 39;

 end
     set(f,'CurrentAxes',ax(end));

 h_xlab = xlabel('Year');
 h_xlab.Units = 'Normalized';
 h_xlab.Position(2) = 0;
 set(gcf, 'Units','centimeter','outerposition', [.25 .25 [25 32]-0.5]);
set(gcf,'papersize',[25 32 ]);
 orient('landscape')
print(f, sprintf('%s/S9_densification_snow_melt',OutputFolder), '-dtiff')
print(f, sprintf('%s/S9_densification_snow_melt',OutputFolder), '-dpdf','-r0')

%% Statistical analysis of densification
f = figure;
hold on

 for ii =1:length(folderlist)
    load(strcat(folderlist{ii},'/run_param.mat'))
    c.OutputFolder = folderlist{ii};
    % extract surface variables
    namefile = sprintf('%s/track_density.nc',folderlist{ii});
    finfo = ncinfo(namefile);
    names={finfo.Variables.Name};
    for i= 1:length(names)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end
    
      % extract run parameters
    load(strcat(folderlist{ii},'/run_param.mat'))
    c.OutputFolder = folderlist{ii};
    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',folderlist{ii},1);
    %             finfo = ncinfo(namefile);
    %             names={finfo.Variables.Name};
    names = {'Year' 'Day' 'Hour' 'H_surf' 'SRout_mdl' 'LRout_mdl' 'SHF' 'LHF' ...
         'GF' 'rainHF' 'meltflux'   'runoff' 'snowfall' 'rainfall'...
          'sublimation' 'snowbkt'};
    for i= 1:length(names)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end

    time_mod = datenum(Year,1,Day,Hour,0,0);
    
     [~ , year, day, hour, pres,...
    T1, T2, z_T1, z_T2, o_T1,o_T2, ...
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, ...
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_out, c]...
    = ExtractAWSData(c);


    runoffhour = -[runoff(1); runoff(2:end)-runoff(1:end-1)] ;
    SMB = snowfall+rainfall+sublimation+runoffhour;
    time = time_mod;
    SMB_hour = table(time, snowfall, rainfall, sublimation, runoffhour,SMB);
    SMB_hour.Properties.VariableNames{5} = 'runoff';

    SEB_hour = table(time_mod,SHF,LHF,SRin - SRout,LRin - LRout_mdl,...
        rainHF,GF, meltflux,meltflux.*c.dt_obs./c.dev./c.L_fus./c.rho_water,...
        'VariableNames',{'time','SHF','LHF','SRnet','LRnet','RainHF','GF','MeltEnergy',...
        'Melt_mweq'});
    
    rho_avg = density_avg_20(6,:);
    delta_rho_tot = 0*rho_avg;
    delta_rho_tot(2:end) = rho_avg(2:end) - rho_avg(1:end-1);
    pore_space = 20 - rho_avg * 20 / c.rho_ice;

    delta_rho_comp = 0*rho_avg;
    delta_rho_comp(2:end) = density_avg_20(1,2:end) - rho_avg(1:end-1);
    delta_rho_precip =  density_avg_20(2,:) -  density_avg_20(1,:);
    delta_rho_subl =  density_avg_20(3,:) -  density_avg_20(2,:);
    delta_rho_melt =  density_avg_20(4,:) -  density_avg_20(3,:);
    delta_rho_runoff =  density_avg_20(5,:) -  density_avg_20(4,:);
    delta_rho_rfrz =  density_avg_20(6,:) -  density_avg_20(5,:);
    
    delta_rho_melt = delta_rho_melt + delta_rho_rfrz;
    
    densif_drivers_hour = table();
    densif_drivers_hour.time = time_mod;
    densif_drivers_hour.delta_rho_comp = delta_rho_comp';
    densif_drivers_hour.delta_rho_precip = delta_rho_precip';
    densif_drivers_hour.delta_rho_subl = delta_rho_subl';
    densif_drivers_hour.delta_rho_melt = delta_rho_melt';
    densif_drivers_hour.delta_rho_tot = delta_rho_tot';
    
    varname = densif_drivers_hour.Properties.VariableNames;
    disp(station{ii})
    
    explained_variance = NaN(5,3);
    for i = 2:5
        lm1 = fitlm(densif_drivers_hour.(varname{i}),...
            delta_rho_tot, 'linear');
        explained_variance(i-1,1) =  lm1.Rsquared.Ordinary*100;
    end

    densif_drivers_daily = AvgTable(densif_drivers_hour, 'daily','sum');
    for i = 2:5
        lm1 = fitlm(densif_drivers_daily.(varname{i}),...
            densif_drivers_daily.delta_rho_tot, 'linear');
        explained_variance(i-1,2) =  lm1.Rsquared.Ordinary*100;


    end
    i = 6;
            lm1 = fitlm(densif_drivers_daily.delta_rho_melt+densif_drivers_daily.delta_rho_precip,...
            densif_drivers_daily.delta_rho_tot, 'linear');
        explained_variance(i-1,2) =  lm1.Rsquared.Ordinary*100;


    densif_drivers_yearly = AvgTable(densif_drivers_hour, 'Jun-yearly','sum');
    for i = 2:5
        lm1 = fitlm(densif_drivers_yearly.(varname{i}),...
            densif_drivers_yearly.delta_rho_tot, 'linear');
        explained_variance(i-1,3) =  lm1.Rsquared.Ordinary*100;
                 g=figure('Visible','off');
            lm1.plot	
        title(station{ii})
        ylabel('density change by hour')
        xlabel(varname{i})
        print(g,...
            sprintf('./Output/Standard runs/%s_%s_lm',station{ii},varname{i}),...
            '-dtiff')
        close(g)
    end
    i = 6;
            lm1 = fitlm(densif_drivers_yearly.delta_rho_melt+densif_drivers_yearly.delta_rho_precip,...
            densif_drivers_yearly.delta_rho_tot, 'linear');
        explained_variance(i-1,3) =  lm1.Rsquared.Ordinary*100;
                 g=figure('Visible','off');
            lm1.plot	
        title(station{ii})
        ylabel('density change by hour')
        xlabel('melt + precip contrib')
        print(g,...
            sprintf('./Output/Standard runs/%s_%s_lm',station{ii},'melt+precip'),...
            '-dtiff')
        close(g)
    disp('       Hourly        Daily       Yearly')
    disp(explained_variance)
    
    DV = datevec(densif_drivers_yearly.time);    
    stairs(DV(:,1), densif_drivers_yearly.delta_rho_tot,'LineWidth',2)
    
    disp('----------------------')
            
 end
 
 axis tight
 box on
 set(gca,'XMinorTick','off','YMinorTick','on','XTick',1998:2014)
 legend(station,'Location','NorthWest')
 legend boxoff
 plot([1998 2014],[0 0],'--k')
 
print(f, sprintf('%s/densification_yearly',OutputFolder), '-dtiff')

%% Trend analysis annual values
clc
trend_annual = table;
trend_annual.no = (1:3)';

trend_annual.start_year(1) = 1999;
trend_annual.end_year(1) = 2009;

trend_annual.start_year(2) = 1999;
trend_annual.end_year(2) = 2014;

trend_annual.start_year(3) = 2001;
trend_annual.end_year(3) = 2014;

for ii =1:length(station)
     load(strcat(folderlist{ii},'/run_param.mat'))
    c.OutputFolder = folderlist{ii};

[time, year, day, hour, pres,...
    T1, T2, z_T1, z_T2, o_T1,o_T2, ...
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, ...
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_out, c] = ExtractAWSData(c);
    data_yearly = AvgTableJJA(data_AWS,'mean');
    DV = datevec(data_yearly.time);
    years = DV(:,1);
    T_year = data_yearly.AirTemperatureC;
    
    figure
    scatter(years,T_year)
    hold on
    Plotlm(years,T_year,'Annotation','on')
    title(station{ii})
    name_slope = sprintf('slope_%s',station{ii});
    name_pvalue = sprintf('pvalue_%s',station{ii});
    name_rmse = sprintf('rmse_%s',station{ii});
    if strcmp(station{ii},'Dye-2')
        name_slope = 'slope_DYE2';
        name_pvalue = 'pvalue_DYE2';
        name_rmse = 'rmse_DYE2';
    end
    if strcmp(station{ii},'NASA-SE')
        name_slope = 'slope_NASASE';
        name_pvalue = 'pvalue_NASASE';
        name_rmse = 'rmse_NASASE';
    end
    
    if strcmp(station{ii},'Crawford Point')
        name_slope = 'slope_CP';
        name_pvalue = 'pvalue_CP';
        name_rmse = 'rmse_CP';
    end
    
    trend_annual.(name_slope) = (1:length(trend_annual.no))';
    trend_annual.(name_pvalue) = (1:length(trend_annual.no))';
    trend_annual.(name_rmse) = (1:length(trend_annual.no))';
    
    years_nonan = years(~isnan(T_year)); %years for which annual average is available
    
    for j = 1:length(trend_annual.no)
        % we calculate things only if start and end years of the period are
        % available
        if ismember(trend_annual.start_year(j),years_nonan) && ...
                ismember(trend_annual.end_year(j),years_nonan)
            years_period = (trend_annual.start_year(j):trend_annual.end_year(j))';
            T_period = T_year(ismember(years,years_period));
            
            lm = fitlm(years_period, T_period);
            trend_annual.(name_slope)(j) = lm.Coefficients.Estimate(2)*10;
            trend_annual.(name_pvalue)(j) = max(lm.Coefficients.pValue);
            trend_annual.(name_rmse)(j) = lm.RMSE;
        else
            trend_annual.(name_slope)(j) = NaN;
            trend_annual.(name_pvalue)(j) = NaN;
            trend_annual.(name_rmse)(j) = NaN;
        end
    end
end

writetable(trend_annual,'./Output/trend_annual.csv','Delimiter' ,';');

