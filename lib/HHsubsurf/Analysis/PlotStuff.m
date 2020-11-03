function [T_subsurf_mod] = PlotStuff(Tsurf_obs, Tsurf, TT, depth_act,depth_act_save, depth_weq, ...
    T_ice, rho_all, rho, snowc, snic, slwc, rfrz, time_mod, Surface_Height, ...
    snowfall, rainfall, runoff,  H_surf,  H_comp, SHF, LHF, ...
    SRin, SRout, LRin, LRout_mdl, rainHF, meltflux, TT_obs, depth_obs, ...
    T_obs, sublimation, OutputFolder, compaction,...
    thickness_act, thickness_weq, T, GF, dgrain, data_AWS, vis, c)


%% ============= Validation Precipitation =====================
% disp('Plotting precipitation')
%    [~, ~, raw] = xlsread('..\AWS_Processing\Input\GCnet\info\Snow pits\GC-Net_AWS_snow_pits_v2005_bapt.xlsx','Overview','A2:C35');
% % [~, ~, raw] = xlsread('..\AWS\GCnet\Snow pits\GC-Net_AWS_snow_pits_v2005_bapt.xlsx','Overview','A2:C35');
% 
% raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
% cellVectors = raw(:,[1,2]);
% raw = raw(~strcmp(raw(:,3),''),3);
% data = reshape([raw{:}],size(raw));
% pit_data = table;
% pit_data.Station = cellVectors(~strcmp(cellVectors(:,1),''),1);
% pit_data.Date = cellVectors(~strcmp(cellVectors(:,2),''),2);
% pit_data.SWEmmweq = data(:,1);
% clearvars data raw cellVectors;
% 
% ind = strcmp(pit_data.Station,c.station);
% pit_data = pit_data(ind,:);
% 
% date_pit = datenum(pit_data.Date);
% ind = and(date_pit>time_mod(1),date_pit<time_mod(end));
% pit_data = pit_data(ind,:);
% 
% if size(pit_data,1)>0
%     col = linspecer(size(pit_data,1),'qualitative');
%     f = figure('Visible','on'); %,'Position',[-2000 100 1800 700]);
%     ha = tight_subplot(1,2,0.05,0.05,0.07);
%     set(f,'CurrentAxes',ha(1))
%     hold on
%     set(f,'CurrentAxes',ha(2))
%     hold on
%     plot(time_mod,Surface_Height,'LineWidth',2)
%     SWE = [];
%     for i=1:size(pit_data,1)
%         date_end = datenum(pit_data.Date(i));
%         [~, ind_end] = min(abs(date_end - time_mod));
%         DV = datevec(date_end);
%         date_start = datenum(DV(:,1)-1,09,01);
%         [~, ind_start] = min(abs(time_mod-date_start));
%         date_start = time_mod(ind_start);
%         
%         SWE(i) = sum(snowfall(ind_start:ind_end))*1000 ...
%             + sum(sublimation(ind_start:ind_end))*1000;
%         set(f,'CurrentAxes',ha(1))
%         scatter(SWE(i),pit_data.SWEmmweq(i),90,'filled','MarkerFaceColor',col(i,:))
%         
%         set(f,'CurrentAxes',ha(2))
%         scatter([date_start date_end],...
%             [Surface_Height(ind_start) Surface_Height(ind_end)],90,'filled','MarkerFaceColor',col(i,:))
%         fprintf('%f \t %f\n',SWE(i),pit_data.SWEmmweq(i))
%         temp =pit_data.Date{i};
%         leg_text{i} = temp((length(temp)-3):length(temp));
%     end
%     bias = nanmean(SWE'-pit_data.SWEmmweq);
%     RMSE = sqrt(nanmean((SWE'-pit_data.SWEmmweq).^2));
%     set(f,'CurrentAxes',ha(1))
%     plot([0 400],[0 400],'k')
%     axis tight square
%     box on
%     set(gca,'XMinorTick','on','YMinorTick','on');
%     legend(leg_text,'Location','NorthWest')
%     xlabel('SWE from station (mm weq)')
%     ylabel('SWE from snow pit (mm weq)')
%     title(sprintf('no correction \nME: %0.2f (%0.1f %%) RMSE: %0.2f', bias, bias/mean(SWE)*100,RMSE))
%     
%     set(f,'CurrentAxes',ha(2))
%     title(c.station)
%     box on
%     datetick('x','mm-yyyy')
%     axis tight square
%     set(gca,'XMinorTick','on','YMinorTick','on','YAxisLocation','right');
%     xlabel('Date')
%     ylabel('Surface Height (m)')
%     
%     if (c.verbose==0)
%         close(f);
%     else
%         print(f, sprintf('%s/precipitation_val',OutputFolder), '-djpeg')
%     end
% end

%% ============ Melt Season Study ==============================
disp('Plotting melt intensity')
data_temp =table(time_mod,meltflux);
data_temp.Properties.VariableNames = {'time','meltflux'};

data_temp = AvgTable(data_temp,'daily',@sum);
melt_daily = data_temp.meltflux;
time_daily = data_temp.time;
melt_season = table();

melting = melt_daily>1;

DV =datevec(data_temp.time);
years_uni = unique(DV(:,1));
melt_info =table();
for i = 1:length(years_uni)
    ind_year = (DV(:,1) == years_uni(i));
    melt_year{i} = melt_daily(ind_year);
    time_year{i}  = data_temp.time(ind_year);
    melt_info.year(i) = years_uni(i);
    melt_info.date_start_year(i) = datenum(years_uni(i),1,1);
    if sum(melt_year{i})>0
        melt_info.start_date(i) = time_year{i}(find(melt_year{i},1,'first'));
        melt_info.end_date(i) = time_year{i}(find(melt_year{i},1,'last'));
        melt_info.duration(i) = melt_info.end_date(i)-melt_info.start_date(i);
        melt_info.max_melt(i) = max( melt_daily(ind_year));
        melt_info.cum_melt(i) = sum( melt_daily(ind_year));
    else
        melt_info.start_date(i) = time_year{i}(10);
        melt_info.end_date(i) = time_year{i}(10);
        melt_info.duration(i) = NaN;
        melt_info.max_melt(i) =NaN;
        melt_info.cum_melt(i) = NaN;
    end
end

if strcmp(c.station,'NUK_K')
    DIV = 10;
else
    DIV = 2;
end

col = lines(length(years_uni));
set(0,'defaultfigurepapersize',[16 29.7]);
set(0,'defaultfigurepaperposition',[.25 .25 fliplr([29.7 16])-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 fliplr([29.7 16])-0.5]);
f = figure('Visible',vis);
hold on
for i = 1:length(years_uni)
    plot([melt_info.start_date(i); melt_info.end_date(i)]'...
        -[melt_info.date_start_year(i)' melt_info.date_start_year(i)'],...
        [melt_info.year(i); melt_info.year(i)],'-o','LineWidth',1,'Color',col(i,:))
    ind = and(time_year{i}>=melt_info.start_date(i),time_year{i}<=melt_info.end_date(i));
    plot(time_year{i}(ind)- melt_info.date_start_year(i), ...
        melt_info.year(i) -  melt_year{i}(ind)/1000/DIV,'LineWidth',2,'Color',col(i,:))
    
    
    ah = annotation('arrow',...
        'headStyle','none','HeadLength',10,'HeadWidth',5);
    set(ah,'parent',gca)
    set(ah,'position',...
        [time_year{i}(find(ind,1,'first'))- melt_info.date_start_year(i),...
        melt_info.year(i),0,-0.928/2],'Color','k');
    plot(time_year{i}(find(ind,1,'first'))- melt_info.date_start_year(i),...
        melt_info.year(i)-0.928/2,'^k','MarkerFaceColor','k')
%     if i==length(years_uni)
    ah2 = annotation('textbox',[1 1 1 1],...
        'String',sprintf('Melt \n%i mm',10*DIV/2),...
        'HorizontalAlignment','Center', 'LineStyle','none');
    set(ah2,'parent',gca)
    set(ah2,'position',...
        [time_year{i}(find(ind,1,'first'))-30 - melt_info.date_start_year(i),...
        melt_info.year(i)-1.5+0.1 1 1],...
        'FitBoxToText','on');
%     end
end
clearvars lm
lm{1} = fitlm(melt_info.year,melt_info.start_date - melt_info.date_start_year);
lm{2} = fitlm(melt_info.year,melt_info.end_date - melt_info.date_start_year);
% lm{3} = fitlm(melt_info.year,melt_info.end_date - melt_info.date_start_year);

% for j = 1:2
%     max(lm{j}.Coefficients.pValue)
%     plot(melt_info.year*lm{j}.Coefficients.Estimate(2) +lm{j}.Coefficients.Estimate(1),...
%         melt_info.year,':k')
% end
axis fill
box on
ylim([years_uni(1)-0.2, years_uni(end)+max( melt_year{i}(ind)/1000/DIV)+0.5])
ylim([years_uni(1)-2, years_uni(end)+0.5])
if strcmp(c.station,'CP1')
    xlim([135 260])
else
    xlim([0 365])
end
set(gca,'XMinorTick','on','YTick',years_uni,'YTickLabel',years_uni,'YDir','reverse')
title(c.station)
xlabel('Day of year')
ylabel('Year')
print(f,sprintf('%s/MeltSeason',OutputFolder),'-djpeg')

set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);

%% ============ Modelled Tsurf vs Obs =========================

if sum(isnan(Tsurf_obs))<length(Tsurf_obs)
    disp('Plotting difference between obs. meas. surface temperature')
    % f = figure('Visible', vis);
    % scatter(Tsurf_obs,Tsurf+c.T_0,'.k')
    % hold on
    % plot(get(gca,'ylim'), get(gca,'ylim'),'k')
    % axis equal
    % xlabel('Tsurf observed')
    % ylabel('Tsurf modelled')
    % print(f, sprintf('%s/Tsurf',OutputFolder), '-djpeg')
    
    temp = datevec(time_mod);
    month_obs = temp(:,2);
    months = ['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; ...
        'Sep'; 'Oct'; 'Nov'; 'Dec'];
    
    g = figure('Visible','off');
    ax = axes('Units','pixels');
    scatter(Tsurf_obs,Tsurf,[],13 - month_obs,'.')
    col =contourcmap('prism',1:12,'colorbar','on','ColorAlignment','center');
    col.YTickLabel = flipud(months);
    hold on
    axis tight
    box on
    plot(get(gca,'ylim'), get(gca,'ylim'),'k')
    set(gca,'XMinorTick','on','YMinorTick','on')
    
    % Create pop-up menu
    popup = uicontrol('Style', 'popup',...
        'String', {'all','Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', ...
        'Sep', 'Oct', 'Nov', 'Dec'},...
        'UserData',struct('Tsurf_obs', Tsurf_obs, 'Tsurf', Tsurf,...
        'month_obs', month_obs, 'c', c),...
        'Position', [20 20 50 20],...
        'Callback', @setmonth);
    
    % Create push button
    btn = uicontrol('Style', 'pushbutton', 'String', 'Clear',...
        'Position', [100 20 50 20],...
        'Callback', 'cla');
    xlabel('Observed surface temperatures (C)')
    ylabel('Modeled surface temperature (C)');
    % Make figure visble after adding all components
    g.Visible = 'on';
    % This code uses dot notation to set properties.
    % Dot notation runs in R2014b and later.
    % For R2014a and earlier: set(f,'Visible','on');
    
    ME = nanmean(Tsurf_obs-Tsurf);
    fprintf('ME = %0.2f\n',ME)
    RMSE = sqrt(nanmean((Tsurf_obs-Tsurf).^2));
    fprintf('RMSE = %0.2f\n',RMSE)
    
    %% box plot month
    f = figure('Visible', 'on');
    ha = tight_subplot(1,2,0.1,0.1,0.1);
    set(f,'CurrentAxes',ha(1));
    
    bias = Tsurf- Tsurf_obs;
    boxplot(bias,month_obs,'Notch','on');
    hold on
    %         axis tight
    ylabel('Surface temperature bias (C)')
    xlimit=get(gca,'XLim');
    set(gca,'XTickLabel',months)
    plot(xlimit,[0 0],':k')
    
    set(f,'CurrentAxes',ha(2));
    day = SRin>5;
    
    boxplot(bias,day,'Notch','on');
    hold on
    %         axis tight
    ylabel('Surface temperature bias (C)')
    xlimit=get(gca,'XLim');
    set(gca,'XTickLabel',{'Night', 'Day'})
    plot(xlimit,[0 0],':k')
    
    print(f, sprintf('%s/Tsurf_validation',OutputFolder), '-djpeg')
    
    %% time series view
    f = figure('Visible',vis);
    ha = tight_subplot(2,1,[.07 .03],[.15 .08],[.1 .08]);
    
    set(f,'CurrentAxes',ha(1))
    
    plot(time_mod, Tsurf)
    hold on
    plot(time_mod,  Tsurf_obs)
    
    legend('Modelled','Observed', 'Location','SouthWest')
    temp = datevec(time_mod)
    
    if abs(time_mod(1)-time_mod(end))>365
        set (gca, 'XTick', [time_mod(1), ...
            datenum(DV(1,1)+1:DV(end,1),1,1)]);
    else
        set (gca, 'XTick', ...
            [time_mod(1), datenum(DV(1,1),temp.Month(1) +1:temp.Month(1)+12,1)]);
    end
    datetick('x','yyyy-mm','keepticks')
    set(gca,'XTickLabel',[],'XMinorTick','on','YMinorTick','on')
    ylabel('Surface temperature (^o C)','Interpreter','tex')
    
    set(f,'CurrentAxes',ha(2))
    plot(time_mod, Tsurf-Tsurf_obs,'Color',[255,165,0]/255)   %RGB('orange'))
    temp = datevec(time_mod);
    
    if abs(time_mod(1)-time_mod(end))>365
        set (gca, 'XTick', [time_mod(1), ...
            datenum(DV(1,1)+1:DV(end,1),1,1)]);
    else
        set (gca, 'XTick', ...
            [time_mod(1), datenum(DV(1,1),temp.Month(1) +1:temp.Month(1)+12,1)]);
    end
    set(gca,'YMinorTick','on','XMinorTick','on')
    
    datetick('x','yyyy-mm','keepticks')
    xlabel('Time')
    ylabel({'Surface temperature', 'difference (^o C)'},'Interpreter','tex')
    print(f, sprintf('%s/Tsurf_diff',OutputFolder), '-djpeg')
    
end

%% =========== Validation with cores =========================
if exist('Core_all.mat', 'file') == 2
    disp('Plotting comparison with density profiles from cores')
    
    load Core_all
    
    switch c.station
        case 'DYE-2'
            NameStation = 'Dye-2';
            %             i_core = [175 130 5 6 31 54 59 53 31];
            i_core = [5 6 31];
        case 'CP1'
            NameStation = 'Crawford Point';
            %             i_core = [33:35 66:74 354 355 544 639];
            i_core = [66:74 639];
        case 'Summit'
            NameStation = 'Summit';
            i_core =  FindCore(Core,'NearestCodeLocation','Summit'); %[41 42 43 44];
            for i = 1:length(i_core)
                disp(i_core(i))
                disp(Core{i_core(i)}.Info.Name)
            end
            i_core([1:4 10:13])= [];
        case 'NASA-SE'
            NameStation = 'NASA-SE';
            i_core = [23 24 25 26];
        case 'GITS'
            NameStation = 'GITS';
            i_core = FindCore(Core,'NearestCodeLocation','Camp Century');

        otherwise
            NameStation = c.station;
            i_core = FindCore(Core,'NearestCodeLocation',NameStation);
    end
    
    
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
        DV = datevec(time_mod(1));
        if Core{i_core(i)}.Info.DateCored.Year < DV(:,1) - 1
            i_remove = [i_remove, i];
        end
    end
    i_core(i_remove) = [];
    num_plot=6;
    f = figure('Visible',vis);
    [ha, ~] = tight_subplot(1, num_plot, 0.01, [0.12 0.01], [0.05 0.01]);
    count = 0;
    switch c.station
        case 'CP1'
            ylim_core = 10;
        case 'DYE-2'
            ylim_core = 15;
        case 'GITS'
            ylim_core = 90;
        otherwise
            ylim_core = 20;
    end
    
    upper_density_meas = [];
    lower_density_meas = [];
    upper_density_mod = [];
    lower_density_mod = [];
    
    for ii = i_core
        count = count+1;
        if count <= num_plot
            set(f,'CurrentAxes',ha(count))
        else
            i_file = 1;
            NameFile = sprintf('%s/CorePlot%i.jpg',OutputFolder, i_file);
            while exist(NameFile, 'file') == 2
                i_file = i_file + 1;
                NameFile = sprintf('%s/CorePlot%i.jpg',OutputFolder,i_file)  ;
            end
            print(f,NameFile,'-djpeg');
            if strcmp(vis,'off')
                close(f);
            end
            f = figure('Visible',vis);
            [ha, ~] = tight_subplot(1, num_plot, 0.01, [0.12 0.01], [0.05 0.01]);
            count = 1;
            set(f,'CurrentAxes',ha(count))
        end
        
        time_core = datenum(Core{ii}.Info.DateCored);
        
        temp = abs(time_mod - time_core);
        [~, ind_time] = min(temp);
        Core_2{1} = Core{ii};
        
        for jk = 2:length(depth_act_save(:,ind_time))
            if thickness_act(jk,ind_time)>0.20
                ind_depth = find(...
                    and(Core_2{1}.Data.Depth < depth_act_save(jk,ind_time)*100,...
                    Core_2{1}.Data.Depth > depth_act_save(jk-1,ind_time)*100));
                Core_2{1}.Data.Density(ind_depth) = nanmean(Core_2{1}.Data.Density(ind_depth));
                if ~isempty(Core_2{1}.Data.Type_perc)
                    Core_2{1}.Data.Type_perc...
                        (min(length(Core_2{1}.Data.Type_perc),ind_depth)) ...
                        = nanmean(Core_2{1}.Data.Type_perc(min(length(Core_2{1}.Data.Type_perc),ind_depth)));
                end
            end
        end
        
        Core_2{2}.Info.Densities = 'y';
        Core_2{2}.Info.DateCored = datetime(datestr(time_mod(ind_time)));
        Core_2{2}.Info.Name = 'Model';
        Core_2{2}.Info.NearestCodeLocation = NameStation;
        
        depth_max = floor(depth_act_save(end,ind_time)*100);
        Core_2{2}.Data.Density = zeros(depth_max,1);
        Core_2{2}.Data.Depth  = [1:depth_max]';
        Core_2{2}.Data.Type = cell(depth_max,1);
        Core_2{2}.Data.Type_perc = zeros(depth_max,1);
        
        for i = length(depth_act_save(:,ind_time)):-1:1
            ind = find(Core_2{2}.Data.Depth<=100*depth_act_save(i,ind_time));
            Core_2{2}.Data.Density(ind) = rho_all(i,ind_time);
            for j = 1:length(ind)
                Core_2{2}.Data.Type{ind(j)} = 'firn';
            end
            Core_2{2}.Data.Type_perc(ind) = snic(i,ind_time)/thickness_weq(i,ind_time)*100;
        end
        
        ind_up = find(and(~isnan(Core{ii}.Data.Density), ...
            Core{ii}.Data.Depth<=500));
        ind_down = find(and(~isnan(Core{ii}.Data.Density),...
            Core{ii}.Data.Depth>500));
        ind = ind_down>length(Core_2{2}.Data.Density);
        ind_down(ind)=[];
        upper_density_meas = [upper_density_meas mean(Core{ii}.Data.Density(ind_up))];
        lower_density_meas = [lower_density_meas mean(Core{ii}.Data.Density(ind_down))];
        
        upper_density_mod = [upper_density_mod mean(Core_2{2}.Data.Density(ind_up))];
        lower_density_mod = [lower_density_mod mean(Core_2{2}.Data.Density(ind_down))];
        
        OverlapPlot(Core_2, [1 2], ...
            'PlotStrat','no',...
            'lag',0,...
            'PlotDifference','no',...
            'YLimit', ylim_core,...
            'span',floor(median(thickness_act(:,ind_time) )*100));
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
    NameFile = sprintf('%s/CorePlot%i.jpg',OutputFolder, i_file);
    while exist(NameFile, 'file') == 2
        i_file = i_file + 1;
        NameFile = sprintf('%s/CorePlot%i.jpg',OutputFolder,i_file)  ;
    end
    print(f,NameFile,'-djpeg');
    
    if strcmp(vis,'off')
        close(f);
    end
    
    f = figure('Visible',vis);
    hold on
    scatter(upper_density_meas, upper_density_mod,120,'filled','MarkerEdgeColor','b','MarkerFaceColor','b')
    scatter(lower_density_meas, lower_density_mod,120,'filled','MarkerEdgeColor','r','MarkerFaceColor','r')
    plot([200 750], [200 750], 'k')
    axis tight square
    box on
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel(sprintf('Measured density (kg/m^3)'))
    ylabel(sprintf('Modelled density (kg/m^3)'))
    legend('Average above 5 m deep','Average below 5 m deep',...
        'Location','SouthWest')
    
    print(f, sprintf('%s/density_comp_2',OutputFolder), '-djpeg')
    if strcmp(vis,'off')
        close(f);
    end
    
    X = table(upper_density_meas', upper_density_mod',...
        lower_density_meas', lower_density_mod');
    X.Properties.VariableNames = {'up_meas','up_mod','down_meas','down_mod'};
    writetable(X,sprintf('%s/%s_density_comp.txt',c.OutputFolder,c.station))
    
    %     if strcmp(c.station, 'CP1')
    %         f = figure('Visible',vis);
    %         [ha, ~] = tight_subplot(1, num_plot, 0.01, [0.12 0.01], [0.05 0.01]);
    %         ii = 206;
    %
    %         set(f,'CurrentAxes',ha(1))
    %         time_core = datenum(Core{ii}.Info.DateCored);
    %         temp = abs(time_mod - time_core);
    %         [~, ind_time] = min(temp);
    %         Core_2{1} = Core{ii};
    %
    %         for jk = 2:length(depth_act_save(:,ind_time))
    %             if thickness_act(jk,ind_time)>0.20
    %                 ind_depth = find(...
    %                     and(Core_2{1}.Data.Depth < depth_act_save(jk,ind_time)*100,...
    %                     Core_2{1}.Data.Depth > depth_act_save(jk-1,ind_time)*100));
    %                 Core_2{1}.Data.Density(ind_depth) = nanmean(Core_2{1}.Data.Density(ind_depth));
    %                 if ~isempty(Core_2{1}.Data.Type_perc)
    %                     Core_2{1}.Data.Type_perc...
    %                         (min(length(Core_2{1}.Data.Type_perc),ind_depth)) ...
    %                         = nanmean(Core_2{1}.Data.Type_perc(min(length(Core_2{1}.Data.Type_perc),ind_depth)));
    %                 end
    %             end
    %         end
    %
    %         Core_2{2}.Info.Densities = 'y';
    %         Core_2{2}.Info.DateCored = datetime(datestr(time_mod(ind_time)));
    %         Core_2{2}.Info.Name = 'Model';
    %         Core_2{2}.Info.NearestCodeLocation = NameStation;
    %
    %         depth_max = floor(depth_act_save(end,ind_time)*100);
    %         Core_2{2}.Data.Density = zeros(depth_max,1);
    %         Core_2{2}.Data.Depth  = [1:depth_max]';
    %         Core_2{2}.Data.Type = cell(depth_max,1);
    %         Core_2{2}.Data.Type_perc = zeros(depth_max,1);
    %
    %         for i = length(depth_act_save(:,ind_time)):-1:1
    %             ind = find(Core_2{2}.Data.Depth<=100*depth_act_save(i,ind_time));
    %             Core_2{2}.Data.Density(ind) = rho_all(i,ind_time);
    %             for j = 1:length(ind)
    %                 Core_2{2}.Data.Type{ind(j)} = 'firn';
    %             end
    %             Core_2{2}.Data.Type_perc(ind) = snic(i,ind_time)/thickness_weq(i,ind_time)*100;
    %         end
    %
    %
    %         OverlapPlot(Core_2, [1 2], ...
    %             'PlotStrat','no',...
    %             'lag',0,...
    %             'PlotDifference','no',...
    %             'YLimit', 50,...
    %             'span',floor(median(thickness_act(:,ind_time) )*100));
    %         xlabel('Density (kg/m^3)')
    %         ylabel('Depth (m)')
    %         for i = 2:num_plot
    %             set(f,'CurrentAxes',ha(i))
    %             set(gca,'Visible','off')
    %         end
    %         i_file = 1;
    %         NameFile = sprintf('%s/CorePlot%i.jpg',OutputFolder, i_file);
    %         while exist(NameFile, 'file') == 2
    %             i_file = i_file + 1;
    %             NameFile = sprintf('%s/CorePlot%i.jpg',OutputFolder,i_file)  ;
    %         end
    %         print(f,NameFile,'-djpeg');
    %         if strcmp(vis,'off')
    %             close(f);
    %         end
    %     end
    
    
    %% Statistical analysis
    
    switch c.station
        case 'DYE-2'
            NameStation = 'Dye-2';
            i_core = [175 130 5 6 31 54 59 53 31];
        case 'CP1'
            NameStation = 'Crawford Point';
            %             i_core = FindCore(Core,'NearestCodeLocation',NameStation);
            i_core = [120 66:74 97 33 199 201 203 205 206 63 52 35];
        otherwise
            NameStation = c.station;
            i_core = FindCore(Core,'NearestCodeLocation',NameStation);
    end
    if ~isempty(i_core)
        for i = 1:length(i_core)
            year_core(i) = Core{i_core(i)}.Info.DateCored.Year;
            depth_range(i) = max(Core{i_core(i)}.Data.Depth);
        end
        if strcmp(c.station,'CP1')
            depth_avg = [300 400 950];
        else
            depth_avg = [300 500 1000];
        end
        year_uni = unique(year_core);
        density_avg = NaN(length(depth_avg),length(year_uni));
        density_std = NaN(length(depth_avg),length(year_uni));
        density_count = NaN(length(depth_avg),length(year_uni));
        for i = 1:length(year_uni)
            ind_year = (year_uni(i) == year_core);
            for j = 1:length(depth_avg)
                ind_core =  (depth_range >= depth_avg(j));
                mean_dens = [];
                ind_combined = find(and(ind_core,ind_year));
                for k = 1:length(ind_combined)
                    mean_dens(k) = nanmean(...
                        Core{i_core(ind_combined(k))}.Data.Density(100:depth_avg(j)));
                end
                density_avg(j,i) = nanmean(mean_dens);
                density_std(j,i) = nanstd(mean_dens);
                density_count(j,i) = sum(~isnan(mean_dens));
            end
        end
        
        f = figure('Visible',vis);
        hold on
        for i = 1:size(density_avg,1)
            errorbar(year_uni(~isnan(density_avg(i,:)))+ (i-1)*0.2,...
                density_avg(i,~isnan(density_avg(i,:))),...
                density_std(i,~isnan(density_avg(i,:))),...
                '-o','LineWidth',2)
        end
        xlabel('Year')
        ylabel('Average density (kg/m^3)')
        for i = 1:length(depth_avg)
            leg_text{i} = sprintf('over 1 - %0.1f m depth range',depth_avg(i)/100);
        end
        legend(leg_text,'Location','north')
        axis tight square
        box on
        set(gca,'XMinorTick','on','YMinorTick','on')
        print(f,sprintf('%s/EvolutionDensity',OutputFolder),'-djpeg');
    end
end

%% =========== Depth of the station ============================
disp('Plotting depth evolution of station mast')
depth_from_surf = depth_act_save;
[~, ind_bot] = min(abs( 4 - depth_from_surf(:,1)));

depth_bot_station(1) = depth_act(ind_bot,1);

for i = 1:length(time_mod)-1
    depth_bot_station(i+1) = depth_bot_station(i) + sum(compaction(ind_bot:end, i));
    [~, ind_bot] = min(abs( depth_bot_station(i+1) - depth_act(:,i+1)));
end

f=figure('Visible',vis);
hold on
plot(time_mod, H_surf)
plot(time_mod,-depth_bot_station)

legend('Surface Height', ...
    'Modelled location of the bottom of the station',...
    'Location','West')

time_mod(1) = time_mod(1);
DV = datevec(time_mod);

if abs(time_mod(1)-time_mod(end))>365
    set (gca, 'XTick', [time_mod(1), ...
        datenum(DV(1,1)+1:DV(end,1),1,1)]);
    datetick('x','yyyy','keepticks')
    set_monthly_tick(time_mod);
else
    %     set (gca, 'XTick', ...
    %         [time_mod(1), datenum(DV(1,1),temp.Month(1) +1:temp.Month(1)+12,1)]);
    datetick('x','dd-mm-yyyy','keepticks')
end
set(gca,'XMinorTick','on','YMinorTick','on')

xlabel('Time')
ylabel('Height (m)')
print(f, sprintf('%s/station_depth',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%%  ============== surface height =========================
disp('Plotting surface height')
f = figure('Visible', vis);
ha=tight_subplot(1,1,0.01,[0.25 0.20],0.2);
set(f,'CurrentAxes',ha(1))
[h, l1, l2] = plotyy(time_mod, Surface_Height-Surface_Height(find(~isnan(Surface_Height),1)),...
    time_mod, snowfall);
hold(h(1));
h1 = plot(h(1),time_mod,H_surf-H_surf(1),'Color', 'k', 'LineWidth', 2);
H_surf_comp = H_surf-H_surf(1)+depth_bot_station'-depth_bot_station(1);
h2 = plot(h(1),time_mod,H_surf_comp,'Color', 'r', 'LineWidth', 2);
l1.Color = [0.4 0.3 0.3];
l2.Color = 'b';
l2.LineWidth = 0.3;
l1.LineWidth = 2;
h(2). YColor = l2.Color ;
hold(h(2))
h3 = plot(h(2),time_mod,...
    meltflux*c.dt_obs/c.dev/c.L_fus/c.rho_water,...
    'Color', 'm', 'LineWidth', 0.5);

axis tight
ylim(h(2),[0 max(max(snowfall,meltflux*c.dt_obs/c.dev/c.L_fus/c.rho_water))*1.10])
ylim(h(1),[min(min(H_surf-H_surf(1),Surface_Height-Surface_Height(1)))*1.1 ...
    max(max(H_surf-H_surf(1),Surface_Height-Surface_Height(1)))*1.1])
xlim(h(1),[time_mod(1) time_mod(end)])
xlim(h(2),[time_mod(1) time_mod(end)])
% h(1).YTick = floor(h(1).YLim(1)*10)/10:...
%     0.1:...
%     floor(h(1).YLim(2)*10)/10;
% legendflex(,'intepreter','tex','Location','NorthWest')
legendflex([l1, l2, h3,h1, h2], {'Observed surface height',...
    'Snowfall', ...
    'Melt', ...
    'Modelled surface height (origin at bottom of column)', ...
    'Modelled surface height (origin at bottom of station mast)'},...
    'ref', gcf, ...
    'anchor', {'n','n'}, ...
    'buffer',[0 0], ...
    'ncol',2, ...
    'fontsize',15);

set_monthly_tick(time_mod);
set(gca,'XMinorTick','on','YMinorTick','on','layer','top',...
    'XTickLabelRotation',20)

xlabel('Time')
ylabel(h(1),'Surface height (m)')
ylabel(h(2),'Modelled snowfall and melt (mweq)')
print(f, sprintf('%s/Melt_surfHeight',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% =========== Validation with SR ============================
disp('Plotting surface height for the melting seasons')

if max(melt_daily)>1000
    melt_daily=melt_daily/1000;
    unit_melt='m';
else
    unit_melt = 'mm' ;
end
DV = datevec(time_mod);
years_list = unique (DV(:,1));

f = figure('Visible',vis,'Position',[0 0 30 length(years_list)/20*70]);
set(f, 'PaperPosition', [0 0 30 length(years_list)/20*70])    % can be bigger than screen 
set(f, 'PaperSize', [30 length(years_list)])

ha = tight_subplot(round(length(years_list)/2), ...
    2,[0.05 0.1], [0.1 0.1] , 0.1);

for i= 1: length(years_list)
    
    i1 = find(time_mod==datenum(years_list(i),4,1),1,'first');
    i2 = find(time_mod==datenum(years_list(i),10,1),1,'first');
    if isempty(i2)
        i2 = length(time_mod);
    end
    if isempty(i1)
        i1 = 1;
    end
    if isempty(i2)
        i1 = length(time_mod);
    end
    set(f,'CurrentAxes',ha(i))
    hold on
    [axout,h1,h3]=plotyy(time_mod(i1:i2),H_surf_comp(i1:i2)-H_surf_comp(i1),...
        time_daily,melt_daily);
    ind = i1:i2;
    h2 = plot( time_mod(i1:i2),Surface_Height(ind)-Surface_Height(...
        ind(find(~isnan(Surface_Height(ind)),1,'first'))),...
    'LineWidth',1.5);

    set_monthly_tick(time_mod(i1:i2))
    h3.Color = [96,96,96]/255; h1.Color = 'b'; h2.Color = 'b'; 
    h1.LineStyle = '--';
    h3.LineWidth = 1.5;  h1.LineWidth = 1.5; h2.LineWidth = 1.5;
    axout(1).YAxis.Color =  h1.Color;
    axout(2).YAxis.Color = h3.Color;
    axout(1).YLabel.String='Surface height (m)';
    axout(2).YLabel.String=['Daily melt (',unit_melt,' w.e.)'];
    axout(2).XAxis.Visible = 'off';
    axout(1).XLim =   [datenum(years_list(i),4,1),...
        datenum(years_list(i),10,1)]; axout(2).XLim = axout(1).XLim;
    box on
    axout(1).YLim = [nanmin([Surface_Height(i1:i2)-Surface_Height(i1);...
        H_surf_comp(i1:i2)-H_surf_comp(i1)]),...
        nanmax([Surface_Height(i1:i2)-Surface_Height(i1);...
        H_surf_comp(i1:i2)-H_surf_comp(i1)])];
    axout(2).YLim = [0 nanmax(melt_daily)*1.1];
    
    title(years_list(i));
    
    set(axout(1),'XTick',datenum(datenum(years_list(i),5:9,1)),...
        'XMinorTick','on','YMinorTick','on')
    datetick('x','keepticks','keeplimits')
    if i<length(years_list)-1
        axout(1).XTickLabel='';
    end
end

legendflex([h1 h2 h3],{'modelled height','observed height','modelled melt'},...
    'ncol',1, 'ref',gcf, 'box','off', 'anchor',{'n' 'n'});

print(f, sprintf('%s/validation_SR',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% ========= Comparison of intial/final density and temperature ===========
disp('Plotting densification')
ind = find(depth_act_save(:,1) <= 20);
tot_vol = depth_act_save(ind(end),1);
tot_mass = sum(rho_all(ind,1) .* thickness_act(ind,1));

rho_avg_1 = tot_mass / tot_vol;

ind = find(depth_act_save(:,end) <= 20);
tot_vol = depth_act_save(ind(end),end);
tot_mass = sum(rho_all(ind,end) .* thickness_act(ind,end));

rho_avg_2 = tot_mass / tot_vol;

f = figure('Visible',vis,'outerposition',[1 1 20 15]);
ha = tight_subplot(1,2,0.05,[0.18 0.02],0.1);
set(f,'CurrentAxes',ha(1))
hold on
plot(depth_act_save(:,1),rho_all(:,1),'b','Linewidth',2)
plot([0 20],[rho_avg_1 rho_avg_1],':b','Linewidth',2)
plot(depth_act_save(:,end), rho_all(:,end),'r','Linewidth',2)
plot([0 20],[rho_avg_2 rho_avg_2],':r','Linewidth',2)

% axis tight 
xlim([0 20])
view([90 90])
xlabel('Depth (m)')
ylabel('Density (kg/m^3)')
box on
set(gca,'XMinorTick','on','XTick',0:20,...
    'XTickLabel',{0, '', 2, '', 4, '', 6, '', 8, '', 10, '', 12, '', 14 '',16,'', 18,'',20},...
    'fontsize',18)
legendflex({datestr(time_mod(1),'dd-mmm-yyyy'),sprintf('mean upper 20 m: %0.2f kg/m^3',rho_avg_1),...
    datestr(time_mod(end),'dd-mmm-yyyy'),sprintf('mean upper 20 m: %0.2f kg/m^3',rho_avg_2)},...
    'anchor',{'e','e'},...
    'ref',gcf);


set(f,'CurrentAxes',ha(2))
hold on
plot(depth_act_save(:,1),rho_all(:,1),'b','Linewidth',2)
plot(depth_act_save(:,end), rho_all(:,end),'r','Linewidth',2)
xlim([0 100])
view([90 90])
box on
set(gca,'XMinorTick','on',...
    'fontsize',18,...
    'Units','normalized',...
    'Position', [0.12 0.2 0.1 0.4],...
    'Box' ,'on',...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'FontSize',10);

print(f, sprintf('%s/firn_evolution_1',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% comparison with the core period
switch c.station
    case 'CP1'
        t1 = datenum(1998,5,1);
        t2 = datenum(2007,5,1);
    case 'DYE-2'
        t1 = datenum(1998,5,1);
        t2 = datenum(2013,5,1);
    case 'NASA-SE'
        t1 = datenum(1998,5,1);
        t2 = datenum(2014,12,31);
    case 'Summit'
        t1 = datenum(2001,5,1);
        t2 = datenum(2014,12,31);    
    case 'GITS'
        t1 = datenum(1966,7,1);
        t2 = datenum(2017,7,1);
    otherwise
        t1 = time_mod(1);
        t2 = time_mod(end);
end
[~, ind_start] = min(abs(time_mod-t1));
[~, ind_end] = min(abs(time_mod-t2));

ind = find(depth_act_save(:,ind_start) <= 20);
tot_vol = depth_act_save(ind(end),1);
tot_mass = sum(rho_all(ind,1) .* thickness_act(ind,1));

rho_avg_1 = tot_mass / tot_vol;

ind = find(depth_act_save(:,ind_end) <= 20);
tot_vol = depth_act_save(ind(end),end);
tot_mass = sum(rho_all(ind,end) .* thickness_act(ind,end));

rho_avg_2 = tot_mass / tot_vol;

f = figure('Visible',vis,'outerposition',[1 1 20 15]);
ha = tight_subplot(1,2,0.05,[0.18 0.02],0.1);
set(f,'CurrentAxes',ha(1))
hold on
plot(depth_act_save(:,ind_start),rho_all(:,ind_start),'b','Linewidth',2)
plot([0 20],[rho_avg_1 rho_avg_1],':b','Linewidth',2)
plot(depth_act_save(:,ind_end), rho_all(:,ind_end),'r','Linewidth',2)
plot([0 20],[rho_avg_2 rho_avg_2],':r','Linewidth',2)
% 
% axis tight square
xlim([0 20])
view([90 90])
xlabel('Depth (m)')
ylabel('Density (kg/m^3)')
box on
set(gca,'XMinorTick','on','XTick',0:20,...
    'XTickLabel',{0, '', 2, '', 4, '', 6, '', 8, '', 10, '', 12, '', 14 '',16,'', 18,'',20},...
    'fontsize',18)
legendflex({datestr(time_mod(ind_start),'dd-mmm-yyyy'),sprintf('mean upper 20 m: %0.2f kg/m^3',rho_avg_1),...
    datestr(time_mod(ind_end),'dd-mmm-yyyy'),sprintf('mean upper 20 m: %0.2f kg/m^3',rho_avg_2)},...
    'anchor',{'e','e'},...
    'ref',gcf);


set(f,'CurrentAxes',ha(2))
hold on
plot(depth_act_save(:,1),rho_all(:,1),'b','Linewidth',2)
plot(depth_act_save(:,end), rho_all(:,end),'r','Linewidth',2)
xlim([0 100])
view([90 90])
box on
set(gca,'XMinorTick','on',...
    'fontsize',18,...
    'Units','normalized',...
    'Position', [0.12 0.2 0.1 0.4],...
    'Box' ,'on',...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'FontSize',10);

print(f, sprintf('%s/firn_evolution_2',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end
%% temperature
f = figure('Visible',vis,'outerposition',[1 1 15 15]);
ha = tight_subplot(1,2,0.05,[0.18 0.02],[0.2 0.1]);

set(f,'CurrentAxes',ha(1))
    hold on
    plot(depth_act_save(:,ind_start),T_ice(:,ind_start)-273.15,'b','Linewidth',2)
    % plot([0 20],[rho_avg_1 rho_avg_1],'--r','Linewidth',2)
    plot(depth_act_save(:,ind_end), T_ice(:,ind_end)-273.15,'r','Linewidth',2)
    % plot([0 20],[rho_avg_2 rho_avg_2],'--b','Linewidth',2)

%     axis tight square
    xlim([0 20])
    view([90 90])
    xlabel('Depth (m)')
    ylabel('Firn temperature (^oC)','Interpreter','tex')
    box on
    set(gca,'XMinorTick','on','XTick',0:20,...
        'XTickLabel',{0, '', 2, '', 4, '', 6, '', 8, '', 10, '', 12, '', 14 '',16,'', 18,'',20},...
        'fontsize',18)

legendflex({datestr(time_mod(ind_start),'dd-mmm-yyyy'),...sprintf('mean upper 20 m: %0.2f kg/m^3',rho_avg_1),...
    datestr(time_mod(ind_end),'dd-mmm-yyyy')},...,sprintf('mean upper 20 m: %0.2f kg/m^3',rho_avg_2)},...
    'anchor',{'e','e'},...
    'ref',gcf);

set(f,'CurrentAxes',ha(2))
    hold on
    plot(depth_act_save(:,1),T_ice(:,ind_start)-273.15,'b','Linewidth',2)
    plot(depth_act_save(:,end), T_ice(:,ind_end)-273.15,'r','Linewidth',2)
    xlim([0 100])
    view([90 90])
    box on
    set(ha(2),'XMinorTick','on',...
        'fontsize',18,...
        'Units','normalized',...
        'Position', [0.4 0.2 0.1 0.4],...
        'Box' ,'on',...
        'XAxisLocation','top',...
        'YAxisLocation','left',...
        'FontSize',10);
    ha(2).YAxisLocation = 'left';
    ha(2).XAxisLocation = 'top';

print(f, sprintf('%s/firn_evolution_2',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end


% [P,S] = polyfit(depth_act_save(depth_act_save(:,1)<20,1),rho_all(depth_act_save(:,1)<20,1),2);
% [P_end,S] = polyfit(depth_act_save(depth_act_save(:,end)<20,end),rho_all(depth_act_save(:,end)<20,1),2);
%
% f = figure ('Visible',vis);
% scatter(depth_act_save(:,1),rho_all(:,1))
% hold on
% plot(depth_act_save(:,1),depth_act_save(:,1).^2.*P(1)+depth_act_save(:,1).*P(2)+P(3))
% scatter(depth_act_save(:,end),rho_all(:,end))
% plot(depth_act_save(:,end),...
%     depth_act_save(:,end).^2.*P_end(1)+depth_act_save(:,1).*P_end(2)+P_end(3))
% xlim([0 20])

%% =========== heat flux to infinite subsurface layer ===========
disp('Plotting heatflux through the bottom layer')
deltaT = T_ice(end,:)-c.Tdeep_AWS-c.T_0;
deltaz = thickness_weq(end,:)./2;
k_eff = 0.021 + 2.5e-6*rho_all(end,:).^2 ;

f = figure ('Visible',vis);
[ha, ~] = tight_subplot(2,1, 0.03, 0.05, 0.08);
set(f,'CurrentAxes',ha(1))
plot(time_mod, k_eff.*deltaT./deltaz)
hold on
plot([time_mod(1) time_mod(end)], [0 0], '--k')
set_monthly_tick(time_mod,gca)
title('Heat flux from last layer downward')
ylabel('Instantaneous (W/m^2)')
set(gca,'XTickLabel',[])

set(f,'CurrentAxes',ha(2))
plot(time_mod, cumsum(k_eff.*deltaT./deltaz.*c.dt_obs))
set_monthly_tick(time_mod,gca)
xlabel('Date')
ylabel('Cumulated (J/m^2)')
print(f, sprintf('%s/bot_hf',OutputFolder), '-djpeg')

%% ====================== Compaction 1 ====================
% if exist('Core_all.mat', 'file') == 2
%
%     load DataCompaction
%     sitename=c.station;
%
%     switch c.station
%         case 'CP1'
%             sitename='Crawford';
%     end
%
%     if sum(strcmp(MetadataCompaction.sitename,sitename))>=1
%         data_site = CompactionAnalysis(DataCompaction, MetadataCompaction, sitename);
%         ind_site = find(strcmp(MetadataCompaction.sitename,sitename));
%
%         % Tracking top and bottom of each borehole, mapping there movement
%         f = figure('Visible', 'on');
%         col = PlotTemp(TT,depth_act,compaction.*1000.*24./thickness_act,...
%             'PlotTherm', 'no',...
%             'PlotIsoTherm', 'no',...
%             'ShowLegend','no',...
%             'cmap','parula',...
%             'XLabel','Time',...
%             'YLabel','Depth (m)',...
%             'CLabel','Compaction rate (mm /m of firn/ day)',...
%             'Range',0:0.1:3);
%         plot(time_mod,-H_surf+H_surf(1),'LineWidth',2)
%         ylim([-max(H_surf)-0.5 25])
%         col.TickLength=col.TickLength*5;
%             datetick('x','yyyy','keepticks')
%
%         depth_from_surf = depth_act;
%         for j = 1:size(depth_act,1)
%             depth_from_surf(j,:) = depth_act(j,:) - depth_act(1,:);
%         end
%
%         depth_top_modscale = NaN(length(data_site),length(time_mod));
%         depth_bot_modscale = NaN(length(data_site),length(time_mod));
%         lines_leg =[];
%         items_leg = cell(0);
%         for j = 1:length(data_site)
%             depth_top = abs(MetadataCompaction.FC_borehole_top_from_surface_m(ind_site,j));
%             depth_bot = abs(MetadataCompaction.FC_borehole_bottom_from_surface_m(ind_site,j));
%             length_1 = depth_bot- depth_top;
%             length_check = abs(MetadataCompaction.FC_borehole_initial_length_m(ind_site,j));
%
%             if abs(length_1-length_check)>0.05
%                 disp('Error in compaction dataset.')
%                 disp('-> Inconsistent borehole length.')
%             end
%             temp =num2str(MetadataCompaction.FC_instrument_installation_daynumbers_YYYYMMDD(ind_site,j));
%             date_start = datenum(str2num(temp(1:4)),str2num(temp(5:6)),str2num(temp(7:8)));
%             if date_start > time_mod(end)
%                 continue
%             end
%             fprintf('\nInstrument %i starting on %s.\n\tDepth of top of borehole: %0.2f\n\tLength of borehol: %0.3f.',...
%                 MetadataCompaction.FC_instrument_IDs(ind_site,j),...
%                 MetadataCompaction.FC_instrument_installation_daynumbers_YYYYMMDD(ind_site,j),...
%                 MetadataCompaction.FC_borehole_top_from_surface_m(ind_site,j),...
%                 MetadataCompaction.FC_borehole_initial_length_m(ind_site,j));
%
%             [~, i_time] = min(abs(time_mod-date_start));
%
%             [~, ind_top] = min(abs( depth_top - depth_from_surf(:,i_time)));
%             [~, ind_bot] = min(abs( depth_bot - depth_from_surf(:,i_time)));
%
%             depth_top_modscale(j,i_time) = depth_act(ind_top,i_time);
%             depth_bot_modscale(j,i_time) = depth_act(ind_bot,i_time);
%
%             for i = i_time:length(time_mod)-1
%                 depth_bot_modscale(j, i+1) = depth_bot_modscale(j, i) + sum(compaction(ind_bot:end, i));
%                 [~, ind_bot] = min(abs( depth_bot_modscale(j, i+1) - depth_act(:,i+1)));
%
%                 depth_top_modscale(j, i+1) = max(depth_top_modscale(j, i) + sum(compaction(ind_top:end, i)),-H_surf(i));
%                 [~, ind_top] = min(abs( depth_top_modscale(j, i+1) - depth_act(:,i+1)));
%             end
%
%             h{j} = plot(time_mod,depth_top_modscale(j,:),'Color',RGB(j),'LineWidth',2);
%             plot(time_mod,depth_bot_modscale(j,:),'Color',RGB(j),'LineWidth',2);
%             lines_leg = [lines_leg, h{j}];
%             items_leg{length(items_leg)+1} = sprintf('Instr. #%i',MetadataCompaction.FC_instrument_IDs(ind_site,j));
%
%         end
%         legend(lines_leg, items_leg,'Location','West')
%         title(sitename)
%         print(f, sprintf('%s/compaction_rate',OutputFolder), '-djpeg')
%         if strcmp(vis,'off')
%             close(f);
%         end
%         % Compaction comparison plot
%         for j = 1:length(data_site)
%         distance_compaction_mod = depth_bot_modscale(j,:) - depth_top_modscale(j,:);
%         compaction_rate_mod = [distance_compaction_mod(1:end-1)-distance_compaction_mod(2:end),...
%             distance_compaction_mod(end-1)-distance_compaction_mod(end)];
%
%         temp =num2str(MetadataCompaction.FC_instrument_installation_daynumbers_YYYYMMDD(ind_site,j));
%         date_start = datenum(str2num(temp(1:4)),str2num(temp(5:6)),str2num(temp(7:8)));
%         if date_start > time_mod(end)
%              fprintf('\nInstr %i out of time range.\n',...
%                  MetadataCompaction.FC_instrument_IDs(ind_site,j));
%             continue
%         end
%
%         time_start = max(time_mod(1),data_site{j}.time(1));
%         time_end = min(time_mod(end),data_site{j}.time(end));
%
%         ind_11 = find(time_mod>=time_start,1,'first');
%         ind_21 = find(data_site{j}.time>=time_start,1,'first');
%
%         ind_12 = find(time_mod<=time_end,1,'last');
%         ind_22 = find(data_site{j}.time<=time_end,1,'last');
%
%         f = figure('Visible',vis);
%         ha = tight_subplot(3,1,0.04, [.07 .08], 0.06);
%
%         set(f,'CurrentAxes',ha(1))
%         hold on
%         plot(time_mod(ind_11:ind_12), distance_compaction_mod(ind_11:ind_12),...
%             '--b','LineWidth',2)
%         plot(data_site{j}.time(ind_21:ind_22),data_site{j}.Compaction_Distance_m_2(ind_21:ind_22),...
%             '-b','LineWidth',2)
%
%         plot(time_mod(ind_11:min(ind_12,ind_11+60*24)), distance_compaction_mod(ind_11:min(ind_12,ind_11+60*24)),...
%             '--','Color',RGB('light light blue'),'LineWidth',2)
%         plot(data_site{j}.time(ind_21:min(ind_22,ind_21+60)),data_site{j}.Compaction_Distance_m_2(ind_21:min(ind_22,ind_21+60)),...
%             '-','Color',RGB('light light blue'),'LineWidth',2)
%
%         ylabel('Length of borehole (m)')
%         text_title = sprintf('%s - Instr. #%i - Length: %0.2f m - Depth of top: %0.2f m \n',...
%             sitename , ...
%             MetadataCompaction.FC_instrument_IDs(ind_site,j),...
%             -MetadataCompaction.FC_borehole_initial_length_m(ind_site,j),...
%             -MetadataCompaction.FC_borehole_top_from_surface_m(ind_site,j));
%         title(text_title);
%         axis tight
%         set(gca,'XMinorTick','on','YMinorTick','on','XLim',[time_mod(1),time_mod(end)])
%         set(gca,'XTickLabel',[]);
%         legend('Modelled','Observed','Location','SouthWest')
%
%         set(f,'CurrentAxes',ha(2))
%         hold on
%         plot(time_mod(ind_11:ind_12), compaction_rate_mod(ind_11:ind_12)*1000*24,...
%             '--r','LineWidth',2)
%         plot(data_site{j}.time(ind_21:ind_22),data_site{j}.Compaction_Rate_md_2(ind_21:ind_22)*1000,...
%             'r','LineWidth',2)
%
%         plot(time_mod(ind_11:min(ind_12,ind_11+60*24)), compaction_rate_mod(ind_11:min(ind_12,ind_11+60*24))*1000*24,...
%             '--','Color',RGB('light light red'),'LineWidth',2)
%         plot(data_site{j}.time(ind_21:min(ind_22,ind_21+60)),data_site{j}.Compaction_Rate_md_2(ind_21:min(ind_22,ind_21+60))*1000,...
%             '-','Color',RGB('light light red'),'LineWidth',2)
%         legend('Modelled','Observed','Location','SouthWest')
%         ylabel('Compac. rate (mm/day)')
%         set(gca,'XTickLabel',[]);
%         axis tight
%         set(gca,'XMinorTick','on','YMinorTick','on','XLim',[time_mod(1),time_mod(end)])
%         if  ismember(MetadataCompaction.FC_instrument_IDs(ind_site,j),[4 2])
%             ylim([0 1])
%         end
%
%         set(f,'CurrentAxes',ha(3))
%         slwc_tt = sum(slwc,1);
%         [ax, h1, h2] = plotyy(time_mod(ind_11:ind_12),T(ind_11:ind_12),...
%             time_mod(ind_11:ind_12),slwc_tt(ind_11:ind_12)*1000);
%         h1.LineWidth = 2;
%         h1.Color = 'c';
%         ax(1).YColor = h1.Color;
%         h2.LineWidth = 2;
%         h2.Color = RGB('dark blue');
%         ax(2).YColor = h2.Color;
%
%         ylabel(ax(1),'Air temperature (m)')
%         ylabel(ax(2),'Total LWC (mm_weq)')
%         xlabel('Date')
%         datetick('x','mm-yyyy','keepticks')
%
%         axes(ax(1))
%         axis tight
%         set(gca,'XMinorTick','on','YMinorTick','on','XLim',[time_mod(1),time_mod(end)])
%
%         axes(ax(2))
%         axis tight
%         set(gca,'XMinorTick','on','YMinorTick','on','XLim',[time_mod(1),time_mod(end)])
%
%         i_file = 1;
%         NameFile = sprintf('%s/compaction_comparison_%i',...
%             OutputFolder,i_file);
%         while exist(strcat(NameFile,'.jpg'), 'file') == 2
%             i_file = i_file + 1;
%             NameFile = sprintf('%s_%i',NameFile(1:end-2),i_file);
%         end
%         print(f,NameFile,'-djpeg');
%     end
%     end
% end

%% ============ Modelled Temperature =========================
disp('Plotting modelled temperature ')
range_temp = -30:2:0;
if strcmp(c.station,'NUK_K')||strcmp(c.station,'Miege')
    range_temp = -10:0.5:0;
end
step = 100;
f = figure('Visible', vis);
ha = tight_subplot(2,1,0.02,[0.2 0.02],0.07);
set(ha(1),'Position',[0.1 0.345 0.86 0.6])
set(f,'CurrentAxes',ha(1))
col = PlotTemp(TT(:,1:step:end),depth_act(:,1:step:end),T_ice(:,1:step:end) - c.T_0,...
    'PlotTherm', 'no',...
    'PlotIsoTherm', 'no',...
    'ShowLegend','no',...
    'XLabel','Time',...
    'YLabel','Depth (m)',...
    'CLabel','Modelled Temperature (^oC)',...
    'cmap','jet',...
    'Range', range_temp);
if  strcmp(c.station, 'GITS')
    ylim([-5 30])
else
    plot(time_mod,-H_surf+H_surf(1),'LineWidth',2)
    ylim([min([-H_surf+H_surf(1)])-1 30])
end
col.TickLength=col.TickLength*5;
col.Position = [ .88    0.3300    0.0430    0.6000];
temp = col.YTickLabel;
for i = 1:length(temp)
    if i/2==floor(i/2)
        temp(i,:)=' ';
    end
end
set(gca,'Color',[0.95 0.95 0.95],'XTickLabel','');
xlabel('')

% temperature at 10 m
T_deep_obs = zeros(1,size(depth_act_save,2));
for i = 1:size(depth_act_save,2)
    T_deep_obs(i) = interp1(depth_act_save(:,i), T_ice(:,i) - c.T_0, 10,'linear');
end

set(ha(2),'Position',[0.1 0.17 0.775 0.14])
set(f,'CurrentAxes',ha(2))
hold on
h1=plot(time_mod,T_deep_obs,'LineWidth',2);
lm = fitlm(time_mod,T_deep_obs);
plot(time_mod, ...
    lm.Coefficients.Estimate(1)+lm.Coefficients.Estimate(2)*time_mod,'k')

% annotation(f,'textbox',...
%     [0.11 0.225 0.19 0.071],...
%     'String',{sprintf('Slope: %0.2e deg/a. p-value: %0.2f',...
%     lm.Coefficients.Estimate(2),max(lm.Coefficients.pValue))},...
%     'LineStyle','none','FitBoxToText','on','FontSize',12);
disp('T10m trend')
fprintf('Slope: %0.1f deg/a. p-value: %0.2f\n',...
     lm.Coefficients.Estimate(2)*365, max(lm.Coefficients.pValue))
% plot([time_mod(1) time_mod(end)], [c.Tdeep c.Tdeep]-c.T_0,'--k')
axis fill
% if strcmp(c.station,'NUK_K')
    h2=plot(data_AWS.time, ...
        data_AWS.IceTemperature8C,'r');
    ylim([-8 0])
    xlim(time_mod([1 50]))
% end
legend([h1 h2],'Modelled','Observed','Location','southeast')
xlabel('Year')
ylabel('Temperature at \newline 10 m (deg C)','interpreter','tex')
% set(gca,'XMinorTick','on','YMinorTick','on')
datetick('x','yyyy','keepticks')
xlim([time_mod(1) time_mod(end)])
set_monthly_tick(time_mod,gca);
print(f, sprintf('%s/Temp_mod',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% =================densities =============
disp('Plotting modelled density')
ylimit = max(-H_surf+H_surf(1))*1.1; 

if  time_mod(end)-time_mod(1) <365
    ylimit = 6;
end
if    strcmp(c.station,'Miege')
        ylimit = 40;
end
f = figure('Visible', 'on');
ha = tight_subplot(2,1,0.02,0.07,0.07);
ha(2).Visible='off';
set(f,'CurrentAxes',ha(1))
col = PlotTemp(TT,depth_act,rho_all,...
    'PlotTherm', 'no',...
    'PlotIsoTherm', 'yes',...
    'ValueIsoTherm', c.rho_pco,...
    'ShowLegend','no',...
    'cmap','parula',...
    'XLabel','Time',...
    'YLabel','Depth (m)',...
    'CLabel','Firn density (kg/m^3)',...
    'Range', 300:25:900);
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
col.YTickLabel=temp;

temp = get(gca,'XTickLabel');
for i = 1:length(temp)
    if i/2==floor(i/2)
        temp(i,:)=' ';
    end
end
set(gca,'XTickLabel',temp);

% set(gca,'Color',[0.95 0.95 0.95]);
print(f, sprintf('%s/rho_1',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%%

f = figure('Visible', vis);
col = PlotTemp(TT,depth_act,rho,...
    'PlotTherm', 'no',...
    'PlotIsoTherm', 'no',...
    'ShowLegend','no',...
    'cmap','parula',...
    'XLabel','Time',...
    'YLabel','Depth (m)',...
    'CLabel','Density of snow only (kg/m^3)',...
    'Range', 200:50:900);

if  strcmp(c.station, 'GITS')
    ylim([-5 30])
else
    plot(time_mod,-H_surf+H_surf(1),'LineWidth',2)
    ylim([min([-H_surf+H_surf(1)])-1 max(max(depth_act))])
end

col.TickLength=col.TickLength*5;
temp = col.YTickLabel;
for i = 1:length(temp)
    if i/2==floor(i/2)
        temp(i,:)=' ';
    end
end
set(gca,'Color',[0.95 0.95 0.95]);
print(f, sprintf('%s/rho_3',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% ================== snic ===========================================
disp('Plotting ice content')
f = figure('Visible', vis);
ha = tight_subplot(2,1,0.02,0.07,0.09);
ha(2).Visible='off';
set(f,'CurrentAxes',ha(1))
col = PlotTemp(TT,depth_act,snic./(snic + snowc + slwc)*100,...
    'PlotTherm', 'no',...
    'PlotIsoTherm', 'yes',...
    'DataIsoTherm', rho_all,...
    'ValueIsoTherm', 900,...
    'ShowLegend','no',...
    'cmap','cool',...
    'XLabel','Time',...
    'YLabel','Depth (m)',...
    'CLabel','Ice content (% mass)',...
    'Range', 0:5:50);
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
% set(gca,...'Color',[0.95 0.95 0.95],
%     'YLim',[ylimit(1) 3])
datetick('x','dd-mm-yyyy')
axis tight
% ylim([min([-H_surf+H_surf(1)])-1 5])
legend('','Snow - ice transition','Snow surface','Location','NorthEast')
legend boxoff
set(gca,'layer','top','LineWidth',1.5,'TickLength',[0.012 0.012])
print(f, sprintf('%s/snic',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% =================== rfrz ===================
disp('Plotting refreezing')
f = figure('Visible', vis);
count = 0;
rfrz_daily=NaN(size(rfrz,1),length(1:24:length(time_mod)));
for j=1:24:length(time_mod)
    count = count+1;
    rfrz_daily(:,count) = mean(rfrz(:,j:min(length(time_mod),j+23)),2);
end

    ha = tight_subplot(2,1,0.02,[0 0.2],0.12);
ha(2).Visible='off';
set(f,'CurrentAxes',ha(1))
col = PlotTemp(TT(:,1:24:length(time_mod)),...
    depth_act(:,1:24:length(time_mod)),rfrz_daily*1000,...
    'PlotTherm', 'no',...
    'PlotIsoTherm', 'yes',...
    'Interp','on',...
    'DataIsoTherm', rho_all(:,1:24:length(time_mod)),...
    'ValueIsoTherm', 900,...
    'ShowLegend','no',...
    'cmap','parula',...
    'XLabel','Time',...
    'YLabel','Depth (m)',...
    'CLabel','Instantaneous refreezing (mm w.eq./hr)',...
    'Range', 0:0.01:0.1);
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
ylimit_temp = get(gca,'YLim');
set(gca,...'Color',[0.95 0.95 0.95],
    'YLim',[ylimit_temp(1) 3])
datetick('x','dd-mm-yyyy')
axis tight
ylim([min([-H_surf+H_surf(1)])-1 ylimit])
set(gcf, 'Position', [.25 .25 [29.7 16]-0.5],...
    'PaperSize',[29.7 16]);
set(gca,'layer','top','LineWidth',1.5,'TickLength',[0.025 0.025])
legend('','Snow - ice transition','Snow surface','Location','NorthEast')
legend boxoff
print(f, sprintf('%s/rfrz',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% ================== dgrain ===========================================
disp('Plotting grain size')

f = figure('Visible', vis);
col = PlotTemp(TT,depth_act,dgrain,...
    'PlotTherm', 'no',...
    'PlotIsoTherm', 'no',...
    'ShowLegend','no',...
    'cmap','parula',...
    'XLabel','Time',...
    'YLabel','Depth (m)',...
    'CLabel','Grain size (mm)',...
    'Range', 0.1:0.1:2.4);
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
print(f, sprintf('%s/dgrain',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% ==================== lwc ==================
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
    'Range', 0:0.3:3);
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
% temp = get(gca,'XTickLabel');
% for i = 1:length(temp)
%     if i/2==floor(i/2)
%         temp(i,:)=' ';
%     end
% end
% set(gca,'XTickLabel',temp);
set(gca,'Color',[0.95 0.95 0.95]);
print(f, sprintf('%s/slwc',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% ==================== saturation ==================
disp('Plotting water saturation')
f = figure('Visible', vis);
ha = tight_subplot(2,1,0.07,0.07,0.07);
set(ha(2),'Visible','off')
set(f,'CurrentAxes',ha(1))
col = PlotTemp(TT,depth_act,100*slwc./(snowc * c.rho_water ./ rho).*c.rho_ice./(c.rho_ice - rho),...
    'PlotTherm', 'no',...
    'PlotIsoTherm', 'yes',...
    'ValueIsoTherm',[7 18],...
    'ShowLegend','on',...
    'cmap','cool',...
    'XLabel','Time',...
    'YLabel','Depth (m)',...
    'CLabel','Volumetric water content (%)',...
    'Range', 0:2:20);
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
% temp = get(gca,'XTickLabel');
% for i = 1:length(temp)
%     if i/2==floor(i/2)
%         temp(i,:)=' ';
%     end
% end
% set(gca,'XTickLabel',temp);
set(gca,'Color',[0.95 0.95 0.95]);
print(f, sprintf('%s/saturation',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% =================== snowc ================
disp('Plotting snow content')
f = figure('Visible', vis);
col = PlotTemp(TT,depth_act,snowc*1000,...
    'PlotTherm', 'no',...
    'PlotIsoTherm', 'no',...
    'ShowLegend','no',...
    'cmap','parula',...
    'XLabel','Time',...
    'YLabel','Depth (m)',...
    'CLabel','Snow content (mm weq)',...
    'Range', 0:5:100);
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
print(f, sprintf('%s/snowc',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% ======= Subsurface temperature =================
disp('Plotting observed subsurface temperatures')
if  ~strcmp(c.station, 'GITS')&&~strcmp(c.station, 'MIT')

    T_subsurf_mod = zeros(size(depth_obs));
    %interpolating modelled at observed depth
    depth_therm_2 = depth_obs;
    for i = 1:size(depth_obs,1)
        depth_therm_2(i,:) = depth_therm_2(i,:) - Surface_Height';
        ind_therm_out = find(depth_therm_2(i,:)<0,1,'first');
        T_obs(i,ind_therm_out:end) = NaN;
    %         depth_obs(i,ind_therm_out:end) = NaN;
        %     depth_therm_2(i,ind_therm_out:end) = 0;
    end

    for j = 1:length(time_mod)
        ind = and(~isnan(depth_act_save(:,j)),~isnan(T_ice(:,j)));
        T_subsurf_mod(:,j) = interp1( [0; depth_act_save(ind,j)],...
            [Tsurf(j)+c.T_0; T_ice(ind,j)],  depth_obs(:,j))-c.T_0;
    end

    if c.ConductionModel == 1
        Tsurf_diff = zeros(size(Tsurf));
        T_diff = T_obs-T_subsurf_mod;
        title_text ='Excess heat observed (^o C)';
    else
        T_diff = T_subsurf_mod - T_obs;
        Tsurf_diff = (Tsurf-Tsurf_obs+c.T_0);
        title_text ='Modelled - observed temperature (^o C)';
    end

    f = figure('Visible', vis);
    f.OuterPosition = ([0.9 0.9 [21 29.7]-0.5]);
    ha = tight_subplot(3,1,0.022,[0.15, 0.04],[0.15 0.15]);
    set(f,'CurrentAxes',ha(1))

    col1 = PlotTemp(TT_obs,depth_therm_2,T_subsurf_mod,...
        'PlotTherm', 'yes',...
        'ShowLegend','no',...
        'PlotIsoTherm', 'no',...
        'cmap','jet',...
        'XLabel',' ',...
        'YLabel','Depth (m)',...
        'CLabel','Interpolated modelled temperature (^o C)',...
        'Range', -15:1:0);
    hold on
    plot(time_mod, -Surface_Height,'r','LineWidth',2)
    xlim([time_mod(1), time_mod(end)])
    ylim([min(-Surface_Height)-0.1 ylimit])
    col1.TickLength=col1.TickLength*5;
    for i = 1:length(col1.YTickLabel)
        if i/2==floor(i/2)
            col1.YTickLabel(i,:)=' ';
        end
    end
    ylabel(col1,'Interpolated modelled\newline     temperature (deg C)','Interpreter','tex')
    set(gca,'Color',[0.7 0.7 0.7],'XTickLabel','')
    % Subsurface temperature bias
    T_diff = T_subsurf_mod - T_obs+c.T_0;
    Tsurf_diff = (Tsurf-Tsurf_obs+c.T_0);
    title_text ='Modelled - observed temperature (^o C)';

    %         f = figure('Visible', 'on');
    %     ha = tight_subplot(2,1,0.035,[0.1, 0.02],0.09);
    set(f,'CurrentAxes',ha(2))
    hold off
    col2 = PlotTemp(TT_obs,depth_therm_2,T_obs-c.T_0,...
        'PlotTherm', 'yes',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','jet',...
        'XLabel',' ',...
        'YLabel','Depth (m)',...
        'Range', -15:1:0);
    hold on
    plot(time_mod, -Surface_Height,'r','LineWidth',2)
    xlim([time_mod(1), time_mod(end)])
    ylim([min(-Surface_Height)-0.1 ylimit])
    col2.TickLength=col2.TickLength*5;
    for i = 1:length(col2.YTickLabel)
        if i/2==floor(i/2)
            col2.YTickLabel(i,:)=' ';
        end
    end
    ylabel(col2,'         Observed \newline temperature (deg C)','Interpreter','tex')

    set(gca,'XTickLabel','') %,'Color',[0.95 0.95 0.95])
    set(gca,'Color',[0.7 0.7 0.7])

    set(f,'CurrentAxes',ha(3))
    col3 = PlotTemp(TT_obs,depth_therm_2,T_diff,...
        'PlotTherm', 'yes',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','BWR_cmap',...
        'XLabel','Year',...
        'YLabel','Depth (m)',...
        'CLabel',title_text,...
        'Range', -11:2:11); %zeros(1,length(Tsurf)));
    hold on
    plot(time_mod, -Surface_Height,'r','LineWidth',2)
    ylim([min(-Surface_Height)-0.1 ylimit])
    xlim([time_mod(1), time_mod(end)])
    col3.TickLength=col3.TickLength*5;

    clim = col1.YTick;
    colormap(ha(1),'jet')
    colormap(ha(2),'jet')
    colormap(col1,'jet')
    colormap(col2,'jet')
    col1.Children.CDataMapping='scaled';

    %     colormap(ha(2),BWR_cmap)
    for i = 1:length(col3.YTickLabel)
        if i/2~=floor(i/2)
            col3.YTickLabel(i,:)=' ';
        end
    end
    ylabel(col3,'Modelled - observed \newline temperature (deg C)','Interpreter','tex')
    set_monthly_tick(time_mod);

    set(gca,'Color',[0.7 0.7 0.7],'XTickLabelRotation',20)
    orient(f,'portrait')
    f.PaperSize = [29.7 21];
    f.PaperPosition=f.OuterPosition;

    print(f, sprintf('%s/T_ice_diff',OutputFolder), '-djpeg')
    if strcmp(vis,'off')
        close(f);
    end

    %% ============== observed 10 m firn temperature ==============================
    % if time_mod(end)-time_mod(1) >365
    switch c.station
        case 'Miege'
            depth_deep = 25;
        case 'KAN-U'
            depth_deep = 5;    
        case 'DYE-2_HQ'
            depth_deep = 9;
        otherwise
            depth_deep= 10;
    end


    TT_obs= repmat(time_mod',size(depth_obs,1),1);

    depth_therm_2 = depth_obs;
    for i = 1:size(depth_obs,1)
        depth_therm_2(i,:) = depth_therm_2(i,:) + Surface_Height';
        ind_therm_out = find(depth_therm_2(i,:)<0,1,'first');
        T_obs(i,ind_therm_out:end) = NaN;
%         depth_obs(i,ind_therm_out:end) = NaN;
        %     depth_therm_2(i,ind_therm_out:end) = 0;
    end

    [depth_sorted, ind_sorted] = sort(depth_obs,1);
    T_ice_sorted = T_obs;
    T_deep_obs = zeros(size(time_mod));
    DV = datevec(time_mod);
    if strcmp(c.station,'CP1')
        ind = ismember(DV(:,1),2001:2004);
        depth_sorted(8,ind)=NaN;
    end
    
    for i =1:length(time_mod)
        T_ice_sorted(:,i)=  T_obs(ind_sorted(:,i),i);
        
        ind = and (~isnan(depth_obs(:,i)),~isnan(T_obs(:,i)));
        if sum(ind)>3
            T_deep_obs(i) = interp1(depth_obs(ind,i),T_obs(ind,i),depth_deep,'linear','extrap');
        else
            T_deep_obs(i)=NaN;
        end
    end
    
    T_deep_obs = hampel(T_deep_obs,7*24,0.1);
    
    ind_year = DV(:,1)-DV(1,1)+1;
    year_uni = unique(DV(:,1));
    T_deep_avg = accumarray(ind_year,T_deep_obs,[],@nanmean);
    for i = 1:length(year_uni)
        if sum(isnan(T_deep_obs(ismember(DV(:,1),year_uni(i))))) > 2*30*24
            T_deep_avg (i) = NaN;
        end
        if sum(ismember(DV(:,1),year_uni(i))) < 10*30*24
            T_deep_avg (i) = NaN;
        end
    end
    
    switch c.station
        case {'DYE-2','DYE-2_long','Summit'}
            T_deep_avg(end-3:end)=NaN;
        case 'NASA-SE'
            T_deep_avg(end-5:end)=NaN;
    end
        
    f= figure('Visible',vis);
    ha = tight_subplot(2,1,0.02,[0.1 0.03],[0.09 0.03]);
    set(f,'CurrentAxes',ha(2))
    hold on
    plot(time_mod,T_deep_obs-273,'LineWidth',2)
    
    for i=1:length(T_deep_avg)
        plot([datenum(year_uni(i),1,1) datenum(year_uni(i),12,31)],...
            T_deep_avg([i i])-273,'k','LineWidth',2)
    end
    legend(sprintf('Interpolated at %i m',depth_deep), ...'Deepest thermocouple','Second deepest','Third deepest')
        'Year average','Location','NorthEast')
    axis tight
    switch c.station
        case {'DYE-2','DYE-2_long'}
        ylim([-20 -10 ])
    end
   
    set_monthly_tick(temp);
    set (gca, 'XTickLabelRotation',0);
    xlim(time_mod([1 end]))
    ylabel(sprintf('Firn temperature at\n%0.1f m depth (degC)',depth_deep))
    box on
        xlabel('Year')
    switch c.station
        case 'CP1'
            pos = [0.2 0.5 0.01 0.01];
        otherwise
            pos = [0.2 0.4 0.01 0.01];
    end
    annotation(gcf,'textbox',...
        pos,...
        'String',{sprintf('Mean firn temperature at %0.2f m: %0.2f ^oC',...
        depth_deep,nanmean(T_deep_avg)-273.15)},...
        'linestyle','none',...
        'FontSize',14,...
        'FitBoxToText','on',...
        'FontWeight','bold');

    set(f,'CurrentAxes',ha(1))
    hold on
%     plot(time_mod([1 end]), [0 0],'k','LineWidth',2)
    for i = size(depth_sorted,1):-1:1
        plot(time_mod, depth_sorted(i,:))
    end
    axis tight
    set_monthly_tick(time_mod);
        xlim(time_mod([1 end]))
        set(gca,'Ydir','reverse')
        ylabel('Depth below the\newline      surface (m)','Interpreter','tex')
        lbl = {'T1','T2','T3','T4','T5','T6','T7','T8','T9','T10'};
        legendflex(lbl, 'ref', gcf, ...
            'anchor', {'ne','ne'}, ...
            'buffer',3*[-13 -10], ...
            'ncol',5, ...
            'fontsize',13);
        plot(time_mod([1 end]), [depth_deep depth_deep],'--k','LineWidth',1.5)
        set(gca,'XTickLabel',[])

        print(f, sprintf('%s/T10m_obs',OutputFolder), '-djpeg')
        if strcmp(vis,'off')
            close(f);
        end
        
        disp('The 10m firn temperature was:')
        for i = 1:length(year_uni)
            fprintf('%i %0.2f\n',year_uni(i), T_deep_avg(i)-273);
        end

        
        
    if and(~strcmp(c.station,'NUK_K'),length(year_uni)>2)
        f = figure('Visible',vis);
        lm = fitlm(year_uni,T_deep_avg-273.15);
        plot(lm)
        axis tight square
        xlabel('Year')
        ylabel('10 m temperature (degC)')
        title(sprintf('mean = %0.2f degC \n slope = %0.3f degC/dec R^2 = % 0.2f p_value = %0.2f',...
            nanmean(T_deep_avg-273.15),...
            lm.Coefficients.Estimate(2)*10,...
            lm.Rsquared.Ordinary, max(lm.Coefficients.pValue)))
        print(f, sprintf('%s/T10m_obs_lm',OutputFolder), '-djpeg')
        if strcmp(vis,'off')
            close(f);
        end
    end
end
%%  ================= Compaction =========================
disp('Plotting compaction')
f = figure('Visible', vis);
[ax,~,h2]= plotyy(time_mod(2:end),(H_comp(2:end)-H_comp(1:end-1))*1000*24,...
    time_mod,-H_comp);
h2.LineWidth = 2;
linkaxes(ax,'x')
time_mod(1) = time_mod(1);
DV = datevec(time_mod);

if abs(time_mod(1)-time_mod(end))>365
    set (gca, 'XTick', [time_mod(1), ...
        datenum(DV(1,1)+1:DV(end,1),1,1)]);
    datetick('x','yyyy','keepticks')
    set_monthly_tick(time_mod);
else
    set (gca, 'XTick', ...
        [time_mod(1), datenum(DV(1,1),DV(1,2) +1:DV(1,2)+12,1)]);
    datetick('x','mm-yyyy','keepticks')
end
set(gca,'XMinorTick','on','YMinorTick','on')

xlabel('Time')
ylabel(ax(1),'Compaction rate (mm/day)')
ylabel(ax(2),'Surface lowering due to compaction (m)')
set(ax(1),'YLim',[0 3],'YTick',0:3,'YTickLabel',0:3)
print(f, sprintf('%s/Comp',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%%  =============== SEB 1 =========================
disp('Plotting SEB')

SEB_hour = table(time_mod,SHF,LHF,SRin - SRout,LRin - LRout_mdl,...
    rainHF,GF, meltflux,meltflux.*c.dt_obs./c.dev./c.L_fus./c.rho_water,...
    'VariableNames',{'time','SHF','LHF','SRnet','LRnet','RainHF','GF','MeltEnergy',...
    'Melt_mweq'});

SEB_JJA = AvgTableJJA(SEB_hour,@sum);

if time_mod(end)-time_mod(1) >365
    SEB_wateryear = AvgTable(SEB_hour,'water-yearly');
    SEB_wateryear.time = datestr(SEB_wateryear.time);
    
    writetable(SEB_wateryear, sprintf('%s/SEB_wateryear.txt', OutputFolder))
    
    SEB_year = AvgTable(SEB_hour,'yearly');
    SEB_year.time = datestr(SEB_year.time);
    
    writetable(SEB_year, sprintf('%s/SEB_year.txt', OutputFolder))
    
    f = figure('Visible',vis);
    col = lines(size(SEB_year,2));
    hold on
    for i = 2:size(SEB_year,2)
        h(i-1) = scatter(str2num(SEB_year.time(:,end-3:end)),...
            table2array(SEB_year(:,i)),'MarkerEdgeColor',col(i-1,:),...
            'MarkerFaceColor',col(i-1,:));
        [lm, ~] = Plotlm(str2num(SEB_year.time(:,end-3:end)),...
            table2array(SEB_year(:,i)),'Annotation','off');
    end
    ylabel('Average energy flux (W m^{-2})')
    legendflex(h,{SEB_year.Properties.VariableNames{2:end}},...
     'ref', gcf, ...
   'anchor', {'n','n'}, ...
   'buffer',[0 0], ...
   'nrow',2, ...
   'fontsize',12,...
   'Interpreter','none');
axis tight
box on
set(gca,'layer','top')
print(f, sprintf('%s/SEB_overview',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

end

f = figure('Visible', vis);

DV = datevec(time_mod);
index = DV(:,2)-DV(1,2)+1 +12*(DV(:,1) - DV(1,1));
Y = horzcat(accumarray(index, SRin - SRout) ,...
    accumarray(index, LRin - LRout_mdl), accumarray(index, GF), accumarray(index, SHF),...
    accumarray(index, LHF) ... %smooth(rainHF, 30*24),
    )./1000;
Y1 = NaN(size(Y));
Y1(Y>0) = Y(Y>0);
Y2 = NaN(size(Y));
Y2(Y<=0) = Y(Y<=0);
melt_month = accumarray(index,meltflux)/1000;

colormap('jet')
hold on
h1 = bar(DV(1,2):(DV(1,2)+size(Y,1)-1),...
    Y1,1,'stacked','EdgeColor','none','LineWidth',0.001);

bar(DV(1,2):(DV(1,2)+size(Y,1)-1),...
    Y2,1,'stacked','EdgeColor','none','LineWidth',0.001)

h2 = stairs(DV(1,2)-1:(DV(1,2)+size(Y,1)-2),...
    melt_month,'r','LineWidth',2);

box on
set(gca,'XMinorTick','on','YMinorTick','on',...
    'XTick',1:12:(size(Y,1)-1),...
    'XTickLabel',unique(DV(:,1)))
ax = get(gca);
ax.XAxis.MinorTickValues = 1:(size(Y,1)-1+DV(1,2));

axis tight
xlabel('Time')
ylabel('Monthly average energy input (kW/m^2)')
h_plot= [h1 h2];
legendflex(h_plot,{'Net shortwave radiation',...
    'Net longwave radiation','Conductive heat flux','Sensible heat flux',...
    'Latent heat flux','Energy for melt'},...
     'ref', gcf, ...
   'anchor', {'n','n'}, ...
   'buffer',[0 0], ...
   'nrow',2, ...
   'fontsize',8);
box on
if length(unique(DV(:,1)))>30
    set(gca,'XTickLabelRotation',90);
end
% i1=find(DV(:,1)==2006,1,'first');
% i2=find(DV(:,1)==2013,1,'last');
% xlim([time_mod(i1) time_mod(i2)])
set(gca,'layer','top')
print(f, sprintf('%s/SEB',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%%  Comp with RACMO
file_RACMO = 'C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Code\AWS_Processing\Input\Secondary data\RACMO_3h_AWS_sites.nc';
finfo = ncinfo(file_RACMO);
melt_RACMO = ncread(file_RACMO,'melt');
time_RACMO = ncread(file_RACMO,'time')+datenum(1900,1,1);
melt_RACMO = melt_RACMO(:,31);
melt_RACMO(time_RACMO<SEB_hour.time(1))=[];
time_RACMO(time_RACMO<SEB_hour.time(1))=[];
RACMO = table();
RACMO.time = time_RACMO;
RACMO.melt = melt_RACMO;
RACMO_yr = AvgTable(RACMO,'yearly',@sum);
RACMO_yr = AvgTable(RACMO,'yearly',@sum);

figure
stairs(RACMO_yr.time,RACMO_yr.melt*3,'LineWidth',1.5)
hold on
stairs(datenum(SEB_year.time),SEB_year.Melt_mweq*1000,'LineWidth',1.5)
set_monthly_tick(SEB_hour.time)
ylabel('Hourly melt ( mm w.e)')
axis tight
legend('RACMO2.3p2','Vandecrux et al. 2020','Location','northoutside')

figure
scatter(RACMO_yr.melt*3,SEB_year.Melt_mweq*1000,'fill')
hold on
plot([0 900],[0 900],'k')
box on
grid on
xlabel('RACMO2.3p2')
ylabel('Vandecrux et al., 2020')
title('Simulated annual melt at Dye-2 (mm w.e.)')
axis tight

%% ================ SEB 2 =========================
f = figure('Visible', vis);
Y = horzcat(accumarray(index, GF*3600),  accumarray(index, SHF*3600),...
    accumarray(index, LHF*3600) , accumarray(index, (SRin - SRout)*3600) ,...
    accumarray(index, (LRin - LRout_mdl)*3600),  ...
    accumarray(index, meltflux*3600) ... %smooth(rainHF, 30*24),
    )./1000000;

ylabel_list = {'Conductive heat flux',...
    'Sensible heat', ...
    'Latent heat flux',...
    'Net shortwave radiation' ,...
    'Net longwave radiation',... %'RainHF'
    'Energy available for melt'};

ha = tight_subplot(4,1,0.01, [.1 0.05], 0.09);
% mycmap = brighten(hsv,-0.6);
% temp = 1:10:60;
% col = mycmap(temp(randperm(6)),:);
time_mod(1) = time_mod(1);
DV = datevec(time_mod);
col = linspecer(6);

count = 0;
for i = [4 5 2 3 1 6]
    if ismember(i, [4 2 1 6])
        count =count+1;
        set(f,'CurrentAxes',ha(count))
        hold on
        plot([DV(1,2) (DV(1,2)+size(Y,1)-1)], [0 0],'--k');
        h1 = stairs(DV(1,2):(DV(1,2)+size(Y,1)-1),...
            Y(:,i),'LineWidth',1.5,'Color',col(i,:));
    else
        h2 = stairs(DV(1,2):(DV(1,2)+size(Y,1)-1),...
            Y(:,i),'LineWidth',1.5,'Color',col(i,:));
    end
    
    axis tight
    box on
    %     i1=find(DV(:,1)==2006,1,'first');
    %     i2=find(DV(:,1)==2013,1,'last');
    %     xlim([time_mod(i1) time_mod(i2)])
    if i==5 || i==3
        h = [h1 h2];
        txt = {ylabel_list{i-1},ylabel_list{i}};
    else
        h = h1;
        txt = {ylabel_list{i}};
    end
    if ismember(i, [5 3 1 6])
        legendflex(h,txt, ...
            'ref', gca, ...
            'anchor', {'n','n'}, ...
            'buffer',[0 -4], ...
            'ncol',length(h),...
            'box','off');
    end
    
    ylimit = get(gca,'YLim');
    if i == 6
        ylim(ylimit + [0 2])
    else
        if strcmp(c.station,'DYE-2')
            ylim(ylimit + [0 17])
        else
            ylim(ylimit + [0 10])
        end
    end
    
    if i==6
        set(gca,'XMinorTick','on','YMinorTick','on',...
            'XTick',1:12:(size(Y,1)-1),...
            'XTickLabel',unique(DV(:,1)))
        ax = get(gca);
        ax.XAxis.MinorTickValues = 1:(size(Y,1)-1+DV(1,2));
        xlabel('Time')
    else
        set(gca,'XTickLabel','')
    end
    if i==1
        ylabel('Cumulated monthly energy flux \newline to the surface (MJ/m^2)','Interpreter','tex')
    end
end
print(f, sprintf('%s/SEB_2',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end

%% ================ SMB 1 =======================================
disp('Plotting SMB')
% runoffhour = -[runoff(1); runoff(2:end)-runoff(1:end-1)] ;
runoffhour=-runoff;
SMB = snowfall+rainfall+sublimation_mweq+runoffhour;
time = time_mod;

SMB_hour = table(time, snowfall, rainfall, sublimation_mweq, runoffhour,SMB);
SMB_hour.Properties.VariableNames{5} = 'runoff';
SMB_month = AvgTable(SMB_hour,'monthly');

f = figure('Visible', vis);
ha = tight_subplot(2,1,0.02, [.1 .09], 0.08);
set(f,'CurrentAxes',ha(1));
hold on

for i =2:5
    stairs([SMB_hour.time; SMB_hour.time(end)],...
        [SMB_hour.(SMB_hour.Properties.VariableNames{i})*1000; SMB_hour.(SMB_hour.Properties.VariableNames{i})(end)*1000],...
        'LineWidth',2)
end

axis tight
box on
legendflex({SMB_month.Properties.VariableNames{2:end-1}},...
    'ref', gcf, ...
    'anchor', {'n','n'},...
    'nrow',1)
set_monthly_tick(time_mod);
set(gca,'XTickLabel','')
ylabel('Hourly contribution (mm_weq)')

set(f,'CurrentAxes',ha(2));
plot(SMB_hour.time,cumsum(SMB_hour.SMB)*1000,'k','LineWidth',2)
ylabel('Cumulated SMB (mm w.eq.)')
set_monthly_tick(time_mod)
axis tight
set(gca,'YAxisLocation','right')
print(f, sprintf('%s/SMB2',OutputFolder), '-djpeg')
if strcmp(vis,'off')
    close(f);
end
if time_mod(end)-time_mod(1) >365
    f = figure ( 'Visible',vis);
    
    SMB_year = AvgTable(SMB_hour,'yearly');
    SMB_season = AvgTable(SMB_hour,'seasonaly');
    SMB_wateryear = AvgTable(SMB_hour,'water-yearly');
    
    %     set(f,'CurrentAxes',ha(2))
    stairs([SMB_year.time; SMB_year.time(end)],...
        [SMB_year.SMB;SMB_year.SMB(end)],'k','LineWidth',3);
    hold on
    stairs([SMB_wateryear.time; SMB_wateryear.time(end)],...
        [SMB_wateryear.SMB;SMB_wateryear.SMB(end)],'r','LineWidth',3);
    stairs([SMB_season.time; SMB_season.time(end)],...
        [SMB_season.SMB;SMB_season.SMB(end)],'b','LineWidth',3);
    plot([SMB_hour.time(1) SMB_hour.time(end)],[0 0],':k')
    axis tight
    box on
    
    set (gca, 'XTick', [time_mod(1), ...
        datenum(DV(1,1)+1:DV(end,1),1,1)]);
    datetick('x','yyyy','keepticks')
    set_monthly_tick(SMB_hour.time);
    
    xlim([SMB_hour.time(1) SMB_hour.time(end)]);
    legend('summed on calendar years','summed on water years',...
        'summed on seasons','Location','SouthEast')
    
    ylabel('Mass balance (m_weq)')
    xlabel('Date')
    
    print(f, sprintf('%s/SMB3',OutputFolder), '-djpeg')
    if strcmp(vis,'off')
        close(f);
    end
end

%%  =============== SMB 2 =========================
if time_mod(end)-time_mod(1) >365
    if (runoff(end) > 0.05*length(time_mod)/365)
        some_runoff = 1;
        color = parula(4);
    else
        some_runoff = 0;
        color = parula(3);
    end
    
    f = figure('Visible', vis);
    ha = tight_subplot(3,1,0.02, [.07 .03], 0.08);
    
    % set(f,'CurrentAxes',ha(1));
    % SMB = cumsum(snowfall)+ cumsum(rainfall)+ H_subl -runoff;
    % hold on
    % plot(time_mod, cumsum(snowfall+rainfall),'Color',color(1,:), 'LineWidth',3)
    % hh = area(time_mod,horzcat(SMB, -H_subl, runoff));
    % hh(1).FaceColor = color(2,:);
    % hh(2).FaceColor = color(3,:);
    % if some_runoff
    %     hh(3).FaceColor = color(4,:);
    %     legend('Snowfall', 'Cumulative mass balance', ...
    %         'Sublimation','Runoff', 'Location','NorthWest')
    % else
    %     legend('Snowfall', 'Surface mass balance', ...
    %         'Sublimation', 'Location','NorthWest')
    % end
    % DV = datevec(SMB_hour.time);
    % if abs(time_mod(1)-time_mod(end))>365
    %     set (gca, 'XTick', [time_mod(1), ...
    %         datenum(DV(1,1)+1:DV(end,1),1,1)]);
    %     set_monthly_tick(SMB_hour.time);
    % else
    %     set (gca, 'XTick', ...
    %         [time_mod(1), datenum(DV(1,1),temp.Month(1) +1:temp.Month(1)+12,1)]);
    % end
    % set(gca,'XTickLabel',[],'layer','top')
    % axis tight
    % ylabel('Cumulative surface \newline mass Balance (m weq)','Interpreter','tex')
    
    set(ha(1),'Visible','off')
    set(ha(3),'Visible','off')
    
    set(f,'CurrentAxes',ha(2));
    hold on
    [xd, yd] =  stairs([SMB_wateryear.time; SMB_wateryear.time(end)],...
        [SMB_wateryear.snowfall;SMB_wateryear.snowfall(end)]);
    xd2 = [xd; xd];
    xd2(2:2:end) = xd;
    xd2(1:2:end) = xd;
    yd2 = [yd; yd];
    yd2(2:2:end) = yd;
    yd2(1:2:end) = yd;
    h2 = patch([xd(1); xd2; xd(end)],[0; yd2; 0],color(3,:));
    
    
    [xd, yd]  = stairs([SMB_wateryear.time; SMB_wateryear.time(end)],...
        [SMB_wateryear.SMB;SMB_wateryear.SMB(end)]);
    xd2 = [xd; xd];
    xd2(2:2:end) = xd;
    xd2(1:2:end) = xd;
    yd2 = [yd; yd];
    yd2(2:2:end) = yd;
    yd2(1:2:end) = yd;
    h3 = patch([xd(1); xd2; xd(end)],[0; yd2; 0],color(2,:));
    
    if some_runoff
        [xd, yd]  = stairs([SMB_wateryear.time; SMB_wateryear.time(end)],...
            [SMB_wateryear.runoff;SMB_wateryear.runoff(end)]);
        xd2 = [xd; xd];
        xd2(2:2:end) = xd;
        xd2(1:2:end) = xd;
        yd2 = [yd; yd];
        yd2(2:2:end) = yd;
        yd2(1:2:end) = yd;
        h4 = patch([xd(1); xd2; xd(end)],[0; yd2; 0],color(4,:));
    end
    
    h1 = stairs([SMB_wateryear.time; SMB_wateryear.time(end)],...
        [SMB_wateryear.snowfall;SMB_wateryear.snowfall(end)]);
    h1.LineWidth =3;
    h1.Color = color(1,:);
    
    h5 =  stairs([SMB_wateryear.time; SMB_wateryear.time(end)],...
        [SEB_wateryear.Melt_mweq; SEB_wateryear.Melt_mweq(end)]);
    h5.LineWidth =3;
    h5.Color = 'r';
    axis tight
    box on
    if abs(time_mod(1)-time_mod(end))>365
        set (gca, 'XTick', [time_mod(1), ...
            datenum(DV(1,1)+1:DV(end,1),1,1)]);
        set_monthly_tick(SMB_hour.time);
    else
        set (gca, 'XTick', ...
            [time_mod(1), datenum(DV(1,1),temp.Month(1) +1:temp.Month(1)+12,1)]);
    end
    title_text = c.station;
    if strcmp(c.station,'CP1')
        title_text = 'Crawford Point';
    end
    title(title_text);
    % xlim([SMB_hour.time(1) SMB_hour.time(end)]);
    % if some_runoff
    %     legend([h(1:2), h(4), h(3)],...
    %         'Snowfall','Sublimation','Runoff','Mass balance',...
    %         'Location','SouthEast')
    % else
    %     legend(h,'Snowfall','Sublimation','Mass balance','Location','SouthEast')
    %
    % end
    
    % set(gca,'XTickLabel',[],'layer','top','YAxisLocation','right')
    set(gca,'layer','top')
    ylabel('        Yearly contributions to \newline surface mass balance(m w.eq.)','Interpreter', 'tex')
    
    if abs(time_mod(1)-time_mod(end))>365
        set (gca, 'XTick', [time_mod(1), ...
            datenum(DV(1,1)+1:DV(end,1),1,1)]);
        datetick('x','yyyy','keepticks')
        set_monthly_tick(SMB_wateryear.time);
    else
        set (gca, 'XTick', ...
            [time_mod(1), datenum(DV(1,1),temp.Month(1) +1:temp.Month(1)+12,1)]);
        datetick('x','mm-yyyy','keepticks')
    end
    axis tight
    xlabel('Year')
    set(gca,'layer','top')
    if some_runoff
        h4.FaceColor = color(4,:);
        legend('Snowfall', 'Cumulative mass balance', ...
            'Sublimation','Runoff', 'Location','NorthWest')
    else
        [leg,~,~,~] = legend([h1 h2 h3 h5],'Snowfall','Sublimation', 'Net surface mass balance','Melt');
        leg.Position = [0.6 0.72 0.2 0.2];
        leg.Box = 'off';
    end
    
    % set(f,'CurrentAxes',ha(3));
    %
    % hold on
    % h(2:4) = area(time_mod,horzcat( sum(slwc,1)',cumsum(sum(rfrz,1)'), runoff));
    % color_2 = winter(3);
    % h(2).FaceColor=color_2(3,:);
    % h(3).FaceColor=color_2(2,:);
    % h(4).FaceColor=color_2(1,:);
    % h(1) = plot(time_mod, cumsum(meltflux*c.dt_obs/c.dev/c.L_fus/c.rho_water),...
    %     'r', 'LineWidth',3);
    % if abs(time_mod(1)-time_mod(end))>365
    %     set (gca, 'XTick', [time_mod(1), ...
    %         datenum(DV(1,1)+1:DV(end,1),1,1)]);
    %     datetick('x','yyyy','keepticks')
    %     set_monthly_tick(temp);
    % else
    %     set (gca, 'XTick', ...
    %         [time_mod(1), datenum(DV(1,1),temp.Month(1) +1:temp.Month(1)+12,1)]);
    %     datetick('x','mm-yyyy','keepticks')
    % end
    % axis tight
    % xlabel('Year')
    % set(gca,'layer','top')
    % ylabel('Liquid water \newline budget (m w.eq.)','Interpreter','tex')
    % if some_runoff
    % legend(h,'Melt','Liquid storage','Refreezing','Runoff', 'Location','NorthWest')
    % else
    % legend(h,'Melt','Liquid storage','Refreezing', 'Location','NorthWest')
    % end
    print(f, sprintf('%s/SMB',OutputFolder), '-djpeg')
    if strcmp(vis,'off')
        close(f);
    end
else
    SMB_wateryear = table;
    SEB_year = table;
end


%% ============= Densification ========================
disp('Plotting densification drivers')
    DensificationStudy(time_mod,depth_act_save,rho_all,...
        compaction,thickness_act,meltflux,sublimation,snowfall,depth_weq,...
        snic, snowc,SEB_year,SMB_wateryear, vis, c)

%% ================ ploting melt ===================
disp('Plotting melt evolution')
if time_mod(end)-time_mod(1) >365
    
    f = figure('Visible',vis);
    hold on
    x = str2num(SEB_year.time(:,end-3:end));
    y = SEB_year.Melt_mweq;
    lm = fitlm(x,y);
    
    scatter(x,y,80,'o','LineWidth',2);
    plot(x, ...
        lm.Coefficients.Estimate(2)*x+lm.Coefficients.Estimate(1),...
        'LineWidth',1.5)
    [~,ci1] = predict(lm,x,...
        'Alpha',0.05,'Simultaneous',true);
    plot(x,ci1,'--r');
    title(sprintf('%s (slope: %0.2f mweq/a pvalue: %0.2f)',...
        c.station,...
        lm.Coefficients.Estimate(2),...
        max(lm.Coefficients.pValue)))
    xlabel('Year')
    ylabel('Melt (m w.eq.)')
    box on
    axis fill square
    set(gca,'XMinorTick','on','YMinorTick','on')
    print(f, sprintf('%s/melt',OutputFolder), '-djpeg')
    if strcmp(vis,'off')
        close(f);
    end
end

%% ==================== validation of accumulation =======================
if time_mod(end)-time_mod(1) >365
    disp('Plotting accumulation and comparison with cores')
    [~, ~, raw] = xlsread('Input\Extra\Accumulation.xlsx','Overview','A2:C11');
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,[2,3]);
    raw = raw(:,1);
    data = reshape([raw{:}],size(raw));
    Accumulation = table;
    Accumulation.sheet = data(:,1);
    Accumulation.station = cellVectors(:,1);
    Accumulation.description = cellVectors(:,2);
    clearvars data raw cellVectors;
    
    ind = find(strcmp(c.station,Accumulation.station));
    
    if ~isempty(ind)
        accum = {};
        count = 1;
        for i = 1:length(ind)
            name_core{i} = Accumulation.description{ind(i)};
            
            data = xlsread('Input\Extra\Accumulation.xlsx',...
                sprintf('Sheet%i',ind(i)));
            accum{count} = table;
            accum{count}.year = data(:,1);
            accum{count}.SMB = data(:,2);
            clearvars data raw;
            count = count +1;
        end
        
        % change comp_box13 to 1 if you want it to appear in the comparison
        comp_box13 = 0;
        if comp_box13
            % loading Box 2013 accumulation rates
            %             namefile = 'C:\Users\bava\ownCloud\Phd_owncloud\Data\Box 2013\Box_Greenland_Accumulation_annual_1840-1999_ver20140214.nc';
            namefile = '..\Box 2013\Box_Greenland_Accumulation_annual_1840-1999_ver20140214.nc';
            finfo = ncinfo(namefile);
            names={finfo.Variables.Name};
            for i= 1:size(finfo.Variables,2)
                eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
            end
            fprintf('\nData extracted from nc files.\n');
            
            dist = sqrt((lat-c.lat).^2 + (lon-c.lon).^2);
            
            [dist_sorted, ind] = sort((dist(:)));
            accum_site = zeros(160,1);
            for i = 1:4
                [i , j]  = ind2sub(size(dist),ind(i));
                accum_site = accum_site + squeeze(acc(i,j,:));
            end
            accum_site = accum_site/4000;
        end
        
        % change comp_MAR to 1 if you want it to appear in the comparison
        comp_MAR = 0;
        if comp_MAR
            if strcmp(c.station,'CP1')
                filename = '..\RCM\MAR\MARv3.5.2_20CRv2c_CP.txt';
            else
                filename = '..\RCM\MAR\MARv3.5.2_20CRv2c_Dye-2.txt';
            end
            
            delimiter = ' ';
            formatSpec = '%f%f%[^\n\r]';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
            fclose(fileID);
            MAR = table(dataArray{1:end-1}, 'VariableNames', {'year','accum'});
            clearvars filename delimiter formatSpec fileID dataArray ans;
            
            if strcmp(c.station,'CP1')
                MAR(1:77,:) = [];
            elseif strcmp(c.station,'DYE-2')
                MAR(1:75,:) = [];
            end
        end
        
        % prepraring other values
        DV = datevec(SMB_year.time);
        disp('From station')
        disp(mean(SMB_year.SMB))
        
        years_accum = [];
        accum_all = [];
        for i = 1:length(accum)
            years_accum = [years_accum accum{i}.year'];
            accum_all = [accum_all accum{i}.SMB'];
        end
        
        years_uni = unique(years_accum);
        accum_std = [];
        accum_mean = [];
        num_cores = [];
        for i = 1:length(years_uni)
            ind = years_accum==years_uni(i);
            accum_std = [accum_std std(accum_all(ind))];
            accum_mean = [accum_mean nanmean(accum_all(ind))];
            num_cores = [num_cores sum(ind)];
        end
        x = [years_uni; accum_mean; accum_std; num_cores];
        
        accum_mean = table;
        num_max_cores = max(num_cores);
        
        accum_mean.year = x(1,num_cores>=1)'; %==num_max_cores);
        accum_mean.accum = x(2,num_cores>=1)'; %==num_max_cores);
        accum_mean.std = x(3,num_cores>=1)'; %==num_max_cores);
        
        accum_mean(find(accum_mean.year<1984),:) = [];
        
        disp('average')
        disp(nanmean(accum_mean.accum))
        
        % Plotting
        
        %maximum of standard deviation 0.31 and on average 0.11
        % yu = accum{i}.SMB+0.1;
        % yl = accum{i}.SMB-0.1;
        yu = accum_mean.accum + accum_mean.std;
        yl = accum_mean.accum - accum_mean.std;
        
        if strcmp(c.station,'DYE-2')
            yu = accum_mean.accum + 0.11;
            yl = accum_mean.accum - 0.11;
        end
        
        % Plotting
        f = figure('Visible',vis,'units','normalized','outerposition',[0 0 0.8 1]);
        hold on
        
        for i = [1 2 5]
            plot(accum{i}.year,accum{i}.SMB,...
                'LineWidth',2);
        end
        plot(DV(:,1),SMB_year.SMB,'k','LineWidth',2);
        
        
        %%
        %     h1 = fill([accum_mean.year' fliplr(accum_mean.year')], [yu' fliplr(yl')], ...
        %         [.9 .9 .9], 'linestyle', 'none');
        %     h2 = plot(accum_mean.year, accum_mean.accum,...
        %         'r','LineWidth',2);
        %
        %     h5 = plot(accum_mean.year([1 end]),[mean(accum_mean.accum),mean(accum_mean.accum)],...
        %     '--r','LineWidth',1.2);
        %
        %     h3 = plot(DV(:,1),SMB_year.SMB,'k','LineWidth',2);
        %     h4 = plot(DV(:,1)([1 end]),[mean(SMB_year.SMB), mean(SMB_year.SMB)],...
        %         '--k','LineWidth',1.2);
        
        %     if comp_MAR && ~comp_box13
        %         h6 = plot(MAR.year,MAR.accum/1000,'LineWidth',2,'Color','m');
        %         lbl = {sprintf('Mean accumulation from cores',num_max_cores),'mean'...
        %             'Standard deviation in core-derived accumulation',...
        %             'Station-derived accumulation', 'mean','MARv3.5.2\_20CRv2c'};
        %         legend([h2, h5, h1 h3 h4,h6],lbl,'Location','NorthOutside')
        %     end
        %     if comp_box13 && ~comp_MAR
        %         h6 = plot(1976:1999,accum_site(137:end),'LineWidth',2,'Color','c');
        %         lbl = {sprintf('Mean accumulation from cores',num_max_cores),'mean'...
        %             'Standard deviation in core-derived accumulation',...
        %             'Station-derived accumulation', 'mean','MARv3.5.2\_20CRv2c'};
        %         legend([h2, h5, h1 h3 h4,h6],lbl,'Location','NorthOutside')
        %     end
        %     if comp_box13 && comp_MAR
        %         h6 = plot(MAR.year,MAR.accum/1000,'LineWidth',2,'Color','g');
        %         h7 = plot(1976:1999,accum_site(137:end),'LineWidth',2,'Color','c');
        %         lbl = {sprintf('Mean accumulation from cores',num_max_cores),'mean'...
        %             'Standard deviation in core-derived accumulation',...
        %             'Station-derived accumulation', 'mean','MARv3.5.2\_20CRv2c','Box (2013)'};
        %         legend([h2, h5, h1 h3 h4,h6,h7],lbl,'Location','NorthOutside')
        %     end
        %     if ~comp_box13 && ~comp_MAR
        %         lbl = {sprintf('Mean accumulation from cores',num_max_cores),'mean'...
        %             'Standard deviation in core-derived accumulation',...
        %             'Station-derived accumulation', 'mean'};
        %
        %         legend([h2, h5, h1 h3 h4],lbl,'Location','NorthOutside')
        %     end
        % end
        %%
        axis tight
        set(gca,'layer','top')
        box on
        set(gca,'XMinorTick','on','YMinorTick','on')
        xlabel('Year')
        ylabel('Surface Mass Balance (m w.eq.)')
        xlim([1970 2014])
        legend({name_core{[1 2 5]},'derived from station'},'Location','EastOutside')
        print(f, sprintf('%s/accum_2',OutputFolder), '-djpeg')
        
        
        print(f, sprintf('%s/accum',OutputFolder), '-djpeg')
        if strcmp(vis,'off')
            close(f);
        end
    else
        disp('No core available for historical accumulation')
    end
end

%% Accum validation
time_accum =1940:2017;
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

f = figure('Visible',vis,'Outerposition',[1 1 30 18]);
col = hsv(20);
count = 1;
    leg_text{1} = 'Station-derived';

    ind = distance(c.lat,c.lon,metadata.Latitude,metadata.Longitude)<0.5;
    ind=find(ind);
    year={};
    accum={};
    for i = 1:length(ind)
        opts = delimitedTextImportOptions("NumVariables", 2);

        opts.DataLines = [1, Inf];
        opts.Delimiter = ";";

        opts.VariableNames = ["year", "b_mm"];
        opts.VariableTypes = ["double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        accum{i} = readtable(sprintf('../Accumulation/data/%i.csv',ind(i)), opts);
        clear opts
    end
    
    DV = datevec(SMB_year.time);
    years = DV(:,1);
    S_year = SMB_year.SMB;
    
    hold on   
    h(1) = plot(years,S_year,':ok','LineWidth',2,...
        'MarkerFaceColor','k','MarkerSize',4);

    for i =1:length(accum)
        disp(metadata.Name(ind(i)))

        leg_text{count+1} = metadata.Name{ind(i)};
        [year_ordered, ind_o] = sort(accum{i}.year);
        h(count+1) = plot(year_ordered, accum{i}.b_mm(ind_o),...
            ':o','LineWidth',2,'Color',col(count,:),...
            'MarkerFaceColor',col(count,:),'MarkerSize',4);
        count = count+1;
%         pause
    end  
 
    xlim(time_accum([1 end]))


    set(gca,'layer','top')

   set (gca,'YMinorTick','on','XMinorTick','on',...
       'XTick',time_accum(1:5:end),...
       'TickLength', [0.01 0.025].*2)

        set(gca,'XTickLabel',time_accum(1:5:end)','XTickLabelRotation',45)
       xlabel('Year')

       h_label = ylabel( sprintf('SMB (m w.eq.)'),'Interpreter','tex');
   box on


unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
ppos(3:4) = pos(3:4);
% pos(1:2) = [1 1];
set(gcf,'paperposition',ppos);
set(gcf,'units',unis);
print(f, sprintf('%s/Accumulation_val',OutputFolder), '-djpeg','-r300')

end
