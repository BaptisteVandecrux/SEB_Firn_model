function PlotStuff_NUK_K(time_mod, Tsurf_obs,Tsurf, T_ice, ...
    depth_act_save, rho_all, Surface_Height,SMB_wateryear, OutputFolder,vis,c)

filename = '.\Side studies\NUK_K\Snow pit summary_bapt.xlsx';
    [~,sheets,~ ] = xlsfinfo(filename);

    date_pit = [];
    f = figure;
orient(f,'landscape')

    ha = tight_subplot(2,5,[0.03 0.01], [0.12 0.16], [0.075 0.01]); 
for i = 1:length(sheets)
    [~, ~, raw0_0] = xlsread(filename,sheets{i},'A1:F50');
    date_pit(i) = datenum(raw0_0{1,1},'dd-mm-yyyy')+0.5;
    [~, ~, raw0_1] = xlsread(filename,sheets{i},'H1:I50');
    raw = [raw0_0,raw0_1];
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,6);
    raw = raw(:,[1,2,3,4,5,7,8]);

    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells

    data = reshape([raw{:}],size(raw));
    SnowPit{i} = table;
    SnowPit{i}.DepthBotcm = data(:,1);
    SnowPit{i}.DepthTopcm = data(:,2);
    SnowPit{i}.DepthMiddlecm = data(:,3);
    SnowPit{i}.Snowdensitykgm3 = data(:,4);
    SnowPit{i}.SWEmm = data(:,5);
    SnowPit{i}.Icelayers = cellVectors(:,1);
    SnowPit{i}.DepthThemistorcm = data(:,6);
    SnowPit{i}.SnowTemperaturedegC = data(:,7);
    clearvars data raw raw0_0 raw0_1 cellVectors R;
    
    SnowPit{i}(1:4,:) = [];
    ind = find(~isnan(SnowPit{i}.DepthBotcm),1,'last')+1;
    SnowPit{i}(ind:end,:) = [];
    
    SnowPit{i}.DepthThemistorcm = -SnowPit{i}.DepthThemistorcm + max(SnowPit{i}.DepthTopcm);
    SnowPit{i}.DepthBotcm = -SnowPit{i}.DepthBotcm + max(SnowPit{i}.DepthTopcm);
    SnowPit{i}.DepthMiddlecm = -SnowPit{i}.DepthMiddlecm + max(SnowPit{i}.DepthTopcm);
    SnowPit{i}.DepthTopcm = -SnowPit{i}.DepthTopcm + max(SnowPit{i}.DepthTopcm);

    [~, i_time] = min(abs(time_mod-date_pit(i)));
    mod_depth = [0; depth_act_save(:,i_time)*100];
    mod_dens= [rho_all(:,i_time); NaN];
    mod_temp = [Tsurf(i_time); T_ice(:,i_time)-273.15];
    
    % plotting
    set(f,'CurrentAxes',ha(i))
    hold on
    stairs(mod_dens,mod_depth,'LineWidth',2)    
    stairs([SnowPit{i}.Snowdensitykgm3; SnowPit{i}.Snowdensitykgm3(end)] , ...
        [SnowPit{i}.DepthTopcm; SnowPit{i}.DepthBotcm(end)] ,'LineWidth',2)
       
    tit_obj = title(datestr(date_pit(i),'dd-mmm-yyyy'));
    tit_obj.Units = 'normalized';
    tit_obj.Position = [ 0.4941    1.3        0];

    if i == 3
        xlabel('Snow density (kg/m^3)')
    end
    if i == 1
        ylab_obj = ylabel('Depth below surface (cm)');
        ylab_obj.Units = 'normalized';
        ylab_obj.Position(2) = ylab_obj.Position(2) - 0.5;       
    else
        set(gca,'YTickLabel','')
    end
    if i == 5
                leg_obj = legendflex({'Modelled','Observed'},'ncol',2,'anchor',{'s' 's'});
        leg_obj.Units = 'normalized';
        leg_obj.Position = [0.7    0.01    0.2401    0.0505];
    end
    axis tight
    box on
    ylim([0 300])
    xlim([200 917])

    set(gca,'XMinorTick','on','YMinorTick','on','ydir','rev','layer','top','XAxisLocation','top','TickLength',[0.06, 0.03])

    set(f,'CurrentAxes',ha(i+5))
    hold on
    plot(mod_temp,mod_depth,'--o','LineWidth',2)
    plot(SnowPit{i}.SnowTemperaturedegC, SnowPit{i}.DepthThemistorcm,'--o','LineWidth',2)
    if i == 3
        xlabel('Snow temperature (degC)')
    end
    if i ~=1
                set(gca,'YTickLabel','')
    end
    axis tight
    xlim([-13 0])
    box on
    ylim([0 300])

    set(gca,'XMinorTick','on','YMinorTick','on','ydir','rev','layer','top','TickLength',[0.06, 0.03])
end
set(gcf,'PaperPositionMode','auto');
print(f, sprintf('%s/comp_snow_pit',OutputFolder), '-dpdf')


%% Snow temperature
[~, ~, raw] = xlsread('.\Side studies\NUK_K\SnowTemperature.xlsx','Ark1','B1:I2946');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
SnowTemperature = reshape([raw{:}],size(raw));
clearvars raw R;

[~, ~, raw, dates] = xlsread('.\Side studies\NUK_K\SnowTemperature.xlsx','Ark1','A1:A2946','',@convertSpreadsheetExcelDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
dates = dates(:,1);

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

time_therm = datetime([dates{:,1}].', 'ConvertFrom', 'Excel');
clearvars raw dates R;

time_therm = datenum(time_therm);
depth_therm = repmat(235 - SnowTemperature(1,:)',1,length(time_mod));
depth_therm=depth_therm/100;
 SnowTemperature(1,:) = [];
 time_therm(1) = [];
 
  SnowTemperature(time_therm<time_mod(1),:) = [];
 time_therm(time_therm<time_mod(1)) = [];
  SnowTemperature(time_therm>time_mod(end),:) = [];
 time_therm(time_therm>time_mod(end)) = [];
 
  SnowTemperature(time_therm==6.9396e+05,:) = [];
  time_therm(time_therm==6.9396e+05,:) = [];

 SnowTemp = NaN(size(depth_therm));
 i_1 = find(time_therm(1) == time_mod);
 i_2 = find(time_therm(end) == time_mod);
SnowTemp(:,i_1:i_2) = SnowTemperature';
ind_isnan = find(isnan(Surface_Height));
ind_isnotnan = find(~isnan(Surface_Height));
Surface_Height(ind_isnan) = interp1(ind_isnotnan ,Surface_Height(ind_isnotnan), ind_isnan);

ind_start=find(time_mod == time_therm(1));
ind_end=find(time_mod == time_therm(end));
depth_therm(:,1:ind_start-1) = NaN;
depth_therm(:,ind_end+1:end) = NaN;
depth_therm_2 = depth_therm;
for i = 1:size(depth_therm,1)
    depth_therm_2(i,:) = depth_therm_2(i,:) + Surface_Height' - Surface_Height(ind_start);
    ind_therm_out = find(depth_therm_2(i,:)<0,1,'first');
    depth_therm(i,ind_therm_out:end) = NaN;
%     depth_therm_2(i,ind_therm_out:end) = 0;
end

%% ======= Interpolated modelled temperature =================
disp('Plotting interpolated modelled subsurface temperatures')
TT_obs= repmat(time_mod',size(depth_therm,1),1);

 T_subsurf_mod = zeros(size(depth_therm));
%interpolating modelled at observed depth
 for j = 1:length(time_mod)
     ind = and(~isnan(depth_act_save(:,j)),~isnan(T_ice(:,j)));
    T_subsurf_mod(:,j) = interp1( [0; depth_act_save(ind,j)],...
        [Tsurf(j)+c.T_0; T_ice(ind,j)],  depth_therm_2(:,j))-c.T_0;
 end
 
    f = figure('Visible', vis);
f.OuterPosition = ([0.9 0.9 [21 29.7]-0.5]);
    ha = tight_subplot(3,1,0.02,[0.07, 0.04],[0.09 0.13]);
    set(f,'CurrentAxes',ha(1))

    col1 = PlotTemp(TT_obs,depth_therm,T_subsurf_mod,...
        'PlotTherm', 'yes',...
        'ShowLegend','no',...
        'PlotIsoTherm', 'no',...
        'cmap','hot',...
        'XLabel',' ',...
        'YLabel','Depth (m)',...
        'CLabel','Interpolated modelled temperature (^o C)',...
        'Range', -15:1:0);
    hold on
    plot(time_mod, -Surface_Height+Surface_Height(ind_start),'r','LineWidth',2)
    xlim([time_therm(1), time_therm(end)])
    ylim([min(-Surface_Height+Surface_Height(ind_start))-0.1 2.5])
    col1.TickLength=col1.TickLength*5;
    for i = 1:length(col1.YTickLabel)
        if i/2==floor(i/2)
            col1.YTickLabel(i,:)=' ';
        end
    end
    ylabel(col1,'Interpolated modelled\newline     temperature (deg C)','Interpreter','tex')
    set(gca,'Color',[0.7 0.7 0.7])
    % Subsurface temperature bias 
        T_diff = T_subsurf_mod - SnowTemp;
        Tsurf_diff = (Tsurf-Tsurf_obs+c.T_0);
        title_text ='Modelled - observed temperature (^o C)';
        
%         f = figure('Visible', 'on');
%     ha = tight_subplot(2,1,0.035,[0.1, 0.02],0.09);
    set(f,'CurrentAxes',ha(2))
        col2 = PlotTemp(TT_obs,depth_therm,SnowTemp,...
        'PlotTherm', 'yes',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','hot',...
        'XLabel',' ',...
        'YLabel','Depth (m)',...
        'Range', -15:1:0);
    hold on
    plot(time_mod, -Surface_Height+Surface_Height(ind_start),'r','LineWidth',2)
    xlim([time_therm(1), time_therm(end)])
    ylim([min(-Surface_Height+Surface_Height(ind_start))-0.1 2.5])
    col2.TickLength=col2.TickLength*5;
    for i = 1:length(col2.YTickLabel)
        if i/2~=floor(i/2)
            col2.YTickLabel(i,:)=' ';
        end
    end
    ylabel(col2,'         Observed \newline temperature (deg C)','Interpreter','tex')
    
    set(gca,'XTickLabel','') %,'Color',[0.95 0.95 0.95])
    set(gca,'Color',[0.7 0.7 0.7])
    
    set(f,'CurrentAxes',ha(3))
    col3 = PlotTemp(TT_obs,depth_therm,T_diff,...
        'PlotTherm', 'yes',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','bwr_cmap',...
        'XLabel','Year',...
        'YLabel','Depth (m)',...
        'CLabel',title_text,...
        'Range', -11:2:11); %zeros(1,length(Tsurf)));
        hold on
    plot(time_mod, -Surface_Height+Surface_Height(ind_start),'r','LineWidth',2)
    ylim([min(-Surface_Height+Surface_Height(ind_start))-0.1 2.5])
    xlim([time_therm(1), time_therm(end)])
    col3.TickLength=col3.TickLength*5;

    clim = col1.YTick;
    colormap(ha(1),'hot')
    colormap(ha(2),'hot')
    colormap(col1,'hot')
    colormap(col2,'hot')
    col1.Children.CDataMapping='scaled';

%     colormap(ha(2),BWR_cmap)
    for i = 1:length(col3.YTickLabel)
        if i/2~=floor(i/2)
            col3.YTickLabel(i,:)=' ';
        end
    end
        ylabel(col3,'Modelled - observed \newline temperature (deg C)','Interpreter','tex')
    set_monthly_tick(time_therm);

    set(gca,'Color',[0.7 0.7 0.7])
    orient(f,'portrait')
    f.PaperSize = [29.7 21];
    f.PaperPosition=f.OuterPosition;
    print(f, sprintf('%s/T_ice_diff_NUK_K',OutputFolder), '-dpng')
if strcmp(vis,'off')
    close(f);
end

%% ======= Comparison to stakes ========================
[~, ~, raw] = xlsread('.\Side studies\NUK_K\Stake mass balance.xlsx','Ark1','A5:Z21');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

data = reshape([raw{:}],size(raw));

% SMB_obs{1} = table;
% SMB_obs{1}.MidElevationm = data(:,1);
% SMB_obs{1}.Aream2 = data(:,2);
% SMB_obs{1}.Winterbalancemweq = data(:,3);
% SMB_obs{1}.Winterbalancem3 = data(:,4);
% SMB_obs{1}.Netbalancemweq = data(:,5);
% SMB_obs{1}.Netbalancem3 = data(:,6);

% SMB_obs{2} = table;
% SMB_obs{2}.MidElevationm = data(:,1);
% SMB_obs{2}.Aream2 = data(:,2);
% SMB_obs{2}.Winterbalancemweq = data(:,7);
% SMB_obs{2}.Winterbalancem3 = data(:,8);
% SMB_obs{2}.Netbalancemweq = data(:,9);
% SMB_obs{2}.Netbalancem3 = data(:,10);

SMB_obs{1} = table;
SMB_obs{1}.MidElevationm = data(:,1);
SMB_obs{1}.Aream2 = data(:,2);
SMB_obs{1}.Winterbalancemweq = data(:,11);
SMB_obs{1}.Winterbalancem3 = data(:,12);
SMB_obs{1}.Netbalancemweq = data(:,13);
SMB_obs{1}.Netbalancem3 = data(:,14);

SMB_obs{2} = table;
SMB_obs{2}.MidElevationm = data(:,1);
SMB_obs{2}.Aream2 = data(:,2);
SMB_obs{2}.Winterbalancemweq = data(:,15);
SMB_obs{2}.Winterbalancem3 = data(:,16);
SMB_obs{2}.Netbalancemweq = data(:,17);
SMB_obs{2}.Netbalancem3 = data(:,18);

SMB_obs{3} = table;
SMB_obs{3}.MidElevationm = data(:,1);
SMB_obs{3}.Aream2 = data(:,2);
SMB_obs{3}.Winterbalancemweq = data(:,19);
SMB_obs{3}.Winterbalancem3 = data(:,20);
SMB_obs{3}.Netbalancemweq = data(:,21);
SMB_obs{3}.Netbalancem3 = data(:,22);

% SMB_obs{6} = table;
% SMB_obs{6}.MidElevationm = data(:,1);
% SMB_obs{6}.Aream2 = data(:,2);
% SMB_obs{6}.Winterbalancemweq = data(:,23);
% SMB_obs{6}.Winterbalancem3 = data(:,24);
% SMB_obs{6}.Netbalancemweq = data(:,25);
% SMB_obs{6}.Netbalancem3 = data(:,26);

clearvars data raw R;


    
    figure
    [ha, ~] = tight_subplot(1, 3, 0.05, 0.05, 0.05);
    year_names = {'2014-15','2015-16','2016-17'};
    col1 = hsv(3);
    col2 = hsv(3);
    col3 = hsv(3);

    for i = 1:2:6;
        axes(ha(1))
        hold on
        h1(round(i/2)) = plot(SMB_obs{round(i/2)}.MidElevationm, ...
            SMB_obs{round(i/2)}.Winterbalancemweq,...
            '--o','LineWidth',2,'Color',col1(round(i/2),:),...
            'MarkerFaceColor',col1(round(i/2),:));
        
        if exist('SMB_wateryear','var')==1
%             SEB_wateryear(end,:) = [];
            years_mod = datetime(datestr(SMB_wateryear.time));
            years_mod = years_mod.Year;
            ind_year = years_mod == str2double([year_names{round(i/2)}(1:4)]);
            if sum(ind_year)>0
                plot(710, SMB_wateryear.snowfall(ind_year),...
                    'x','Color',col1(round(i/2),:),...
                    'LineWidth',2,...
                    'MarkerSize',sqrt(140));
            end
        end
        ylabel('m w.eq.')
        title('Winter balance')
        xlabel('Elevation (m a.s.l.)')
        axis fill square
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on

        axes(ha(2))
        hold on
        h2(round(i/2)) = plot(SMB_obs{round(i/2)}.MidElevationm, SMB_obs{round(i/2)}.Netbalancemweq,...
            '--or','LineWidth',2,'Color',col2(round(i/2),:),...
            'MarkerFaceColor',col2(round(i/2),:));
        plot(SMB_obs{round(i/2)}.MidElevationm,zeros(size(SMB_obs{round(i/2)}.MidElevationm)),'--k');
        axis fill square
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on
        if exist('SMB_wateryear','var')==1
            if sum(ind_year)>0
                plot(710, SMB_wateryear.SMB(ind_year),...
                    'x','Color',col2(round(i/2),:),...
                    'LineWidth',2,...
                    'MarkerSize',sqrt(140));
            end
        end
        ylabel('m w.eq.')
        
        title('Net mass balance')
        xlabel('Elevation (m a.s.l.)')

        axes(ha(3))
        hold on
        h3(round(i/2))=plot(SMB_obs{round(i/2)}.MidElevationm, SMB_obs{round(i/2)}.Netbalancemweq-SMB_obs{round(i/2)}.Winterbalancemweq,...
            '--or','LineWidth',2,'Color',col3(round(i/2),:),...
            'MarkerFaceColor',col3(round(i/2),:));
        axis fill square
        set(gca,'XMinorTick','on','YMinorTick','on')
        box on
        if exist('SEB_wateryear','var')==1
            if sum(ind_year)>0
                plot(710, SMB_wateryear.runoff(ind_year)...
                    +SMB_wateryear.sublimation(ind_year),...
                    'x','Color',col3(round(i/2),:),...
                    'LineWidth',2,...
                    'MarkerSize',sqrt(140));
            end
        end

        ylabel('m w.eq.')
                title('Melt')
        xlabel('Elevation (m a.s.l.)')
    end
    hh1 = plot([NaN NaN],[NaN NaN],'xk','LineWidth',2,'MarkerSize',sqrt(140));

    legend([h3 hh1], ...
        {year_names{:},'Derived from the weather station'},...
        'Location','NorthOutside')
    
   %%
   figure
   plot(years_mod, SMB_wateryear.snowfall,'--o','LineWidth',2)
   hold on
   for i = 1:3
       scatter(years_mod(i), SMB_obs{i}.Winterbalancemweq(1),'fill')
    end
end


 

 
 


