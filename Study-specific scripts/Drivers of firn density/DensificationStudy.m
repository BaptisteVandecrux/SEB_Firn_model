function [h_legend] = DensificationStudy(time_mod, density_avg_20, ...
             vis, c)
% in this part we analyze the densidifcation down to 20m and the
% participation of various processes to that net change

%% Extracting average density and contributions
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
    
%% Pore space analysis

% ind_discretize = discretize(0.01:0.01:20,[0; depth_act_save(:,1)]);
% rho_discretize = rho_all(ind_discretize,1);
% rho_avg_start = mean(rho_discretize)
% 
% ind_discretize = discretize(0.01:0.01:20,[0; depth_act_save(:,end)]);
% rho_discretize = rho_all(ind_discretize,end);
% rho_avg_end = mean(rho_discretize)
% disp('Reduction in pore space using regridding:')
% disp(- (rho_avg_end - rho_avg_start) /( c.rho_pco - rho_avg_start) *100)
% 
% rho_avg = zeros(length(time_mod),1);
% pore_space = zeros(length(time_mod),1);
% 
% for i = 1:length(time_mod)
% %     fprintf('%i\n',floor(i/length(time_mod)*100))
%     ind_discretize = discretize(0.01:0.01:20,[0; depth_act_save(:,i)]);
%     rho_discretize = rho_all(ind_discretize,i);
%     rho_avg(i) = mean(rho_discretize);
%     pore_space(i) =  (c.rho_pco - rho_avg(i))*20/c.rho_ice;
% end

diff_dens_year = table;

temp = datetime(datestr(time_mod));
diff_dens_year.Year = unique(temp.Year);
time_stamp = ...
    datenum(diff_dens_year.Year(1):diff_dens_year.Year(end),9,1);
rho_stamp = NaN(size(time_stamp));
ps_stamp = NaN(size(time_stamp));
for i = 1:length(time_stamp)  
    ind = find(time_mod==time_stamp(i));
    if isempty(ind)
       rho_stamp(i) = NaN;
       ps_stamp(i) = NaN;
    else
       rho_stamp(i) = rho_avg(ind);
       ps_stamp(i) = pore_space(ind);
    end
end
% sum(delta_rho_tot(117077:125860))
if time_mod(end)-time_mod(1) >365

    diff_dens_year.DensityChange = zeros(size(diff_dens_year.Year));
    diff_dens_year.DensityChange_perc = zeros(size(diff_dens_year.Year));
    diff_dens_year.DensityChange(2:end) = rho_stamp(2:end)-rho_stamp(1:end-1);
    tot_densification = sum(diff_dens_year.DensityChange(2:end));
    diff_dens_year.DensityChange_perc(2:end) = ...
        diff_dens_year.DensityChange(2:end)/tot_densification*100;

    [~, ind] = sort(diff_dens_year.DensityChange_perc);
    sorted_dens = (diff_dens_year(flipud(ind),:));
    sorted_dens.CumulativeDensityChange = cumsum(sorted_dens.DensityChange_perc);

    % save(sprintf('Side studies/%s_dens.mat',c.station),'sorted_dens')

%     density_tab = {};
%     name_station = {};
%     f = dir('./Side studies');
%     aux1 = {f.name};
%     aux2 = strfind(aux1,'_dens.mat');
%     for i = 1:length(aux2)
%         if ~isempty(aux2{i})
%             load(sprintf('Side studies/%s',aux1{i}));
%             density_tab{length(density_tab)+1} = sorted_dens;
%             name_station{length(name_station)+1} = aux1{i}(1:aux2{i}-1);
%         end
%     end
% 
%     x = linspace(-15,20);
%     col = linspecer(length(name_station)+1);
% 
%     f = figure('Visible',vis);
%     ha = tight_subplot(1, 2, 0.07, 0.07,0.07);
%     set(f,'CurrentAxes',ha(1))
%     hold on
%     leg_text={};
%     for i = 1:length(name_station)
%         scatter(density_tab{i}.Melt*1000,density_tab{i}.DensityChange, ...
%             120,'filled','MarkerFaceColor',col(i,:))
%         exp_fit = fit(density_tab{i}.DensityChange, density_tab{i}.Melt*1000,'exp1');
%         plot(exp_fit(x),x,'LineWidth',1.5,'Color',col(i,:))
%         leg_text = {leg_text{:}, name_station{i},'best fit'};
%     end
% 
%     xlabel('Melt (mm w.eq.)')
%     ylabel('Density change (kg/m3)')
%     box on
%     axis tight square
%     xlim([-50 1000])
%     legend(leg_text,'Location','SouthEast')
%     set(gca, 'XMinorTick','on','YMinorTick','on');
% 
%     set(f,'CurrentAxes',ha(2))
%     hold on
%     for i = 1:length(name_station)
%         scatter(density_tab{i}.Precip*1000,density_tab{i}.DensityChange, ...
%             120,'filled','MarkerFaceColor',col(i,:))
%         [~, ~] = Plotlm(density_tab{i}.Precip*1000, density_tab{i}.DensityChange,...
%             'Annotation','off',...
%             'Color',col(i,:));
%     end
% 
%     xlabel('Total precipitation (mm w. eq.)')
%     ylabel('Density change (kg/m^3)')
%     box on
%     axis tight square
%     ylim([-15 20])
%     set(gca, 'XMinorTick','on','YMinorTick','on');
%     print(f, sprintf('%s/SensitivityDensification',c.OutputFolder), '-dtiff')
%         close(f);
end

%% putting data in table
time = time_mod;

part_dens_hour = table(time, ...
    delta_rho_comp', ...
    delta_rho_melt', ....
    delta_rho_precip', ...
    delta_rho_subl', ...
    delta_rho_runoff', ...
    delta_rho_tot');
part_dens_hour.Properties.VariableNames = { 'time', 'delta_rho_comp', ...
    'delta_rho_melt', ....
    'delta_rho_precip', ...
    'delta_rho_subl', ...
    'delta_rho_runoff', ...
    'delta_rho_tot'};
part_dens_hour.check = part_dens_hour.delta_rho_comp ...
    + part_dens_hour.delta_rho_melt ...
    + part_dens_hour.delta_rho_precip ...
    + part_dens_hour.delta_rho_subl ...
    + part_dens_hour.delta_rho_runoff;

if time_mod(end)-time_mod(1) >365
    part_dens_year = AvgTable(part_dens_hour,'Jun-yearly', 'sum');
    temp = datetime(datestr(part_dens_year.time));
    part_dens_year.time = temp.Year;
    
    ind = ismember(diff_dens_year.Year,    part_dens_year.time);
    part_dens_year.DensityChange = diff_dens_year.DensityChange(ind);
    part_dens_year.check = part_dens_year.delta_rho_comp ...
        + part_dens_year.delta_rho_melt ...
        + part_dens_year.delta_rho_precip ...
        + part_dens_year.delta_rho_subl ...
        + part_dens_year.delta_rho_runoff;
    
    disp(c.station)
    fprintf('var precip        %0.2f\n',var(part_dens_year.delta_rho_precip))
    fprintf('var comp        %0.2f\n',var(part_dens_year.delta_rho_comp))
    fprintf('var melt        %0.2f\n',var(part_dens_year.delta_rho_melt))
end

%% plot month
if c.plot_month_dens
    time_start = datetime(datestr(time_mod(1)));  
    time_end = datetime(datestr(time_mod(end)));  
    time_new = datenum(time_start);
        i = 1;
        while time_new(end) < datenum(time_end)
           time_new(i) = datenum(time_start.Year, time_start.Month + i-1, 1);
           i = i + 1;
        end
    time_dt = datetime(datestr(time_mod));
    % ind_binned = time_dt.Month-time_dt.Month(1)+1 +12*(time_dt.Year - time_dt.Year(1));
    ind_binned = discretize(time_mod,time_new);
    ind_nan =isnan(ind_binned);

    Y = NaN(max(ind_binned),4);
    count =1;
    for i = [4 1 2 3]
        Y(:,count) = accumarray(ind_binned(~ind_nan), ...
            table2array(part_dens_hour(~ind_nan,i+1)));
        count = count +1;
    end

    Z = [time_new(1:end-1)' Y];
%     dlmwrite(sprintf('%s/monthly_densification_drivers.txt',c.OutputFolder),Z)

    Y1 = NaN(size(Y));
    Y1(Y>0) = Y(Y>0);
    Y2 = NaN(size(Y));
    Y2(Y<=0) = Y(Y<=0);
    delta_rho_sum = accumarray(ind_binned(~ind_nan),part_dens_hour.check(~ind_nan));

    time_plot = time_dt.Month(1):(time_dt.Month(1)+size(Y,1)-1);

    f = figure('Visible', vis);
    ha = tight_subplot(2,1,0.03,0.1,0.1);
    set(f,'CurrentAxes',ha(1))
    colormap('jet')
    hold on
    h1 = bar(time_plot,...
        Y1,1,'stacked','EdgeColor','none','LineWidth',0.001);

    bar(time_plot,...
        Y2,1,'stacked','EdgeColor','none','LineWidth',0.001)

    h2 = stairs(time_plot-0.5,...
        delta_rho_sum,'r','LineWidth',1.5);

    box on
    set(gca,'XMinorTick','on','YMinorTick','on',...
        'XTick',1:12:(size(Y,1)-1),...
        'XTickLabel',unique(time_dt.Year))
    ax = get(gca);
        ax.XAxis.MinorTickValues = 1:(size(Y,1)-1+time_dt.Month(1));

    axis tight
    % xlabel('Time')
    ylabel(sprintf('Contribution to \n density change (kg/m^3)'))
    legendflex([h1 h2],{'Sublimation','Dry compaction', 'Melt',...
        'Precipitation','Total change'}, ...
        'ref', gca,'box','off', ...
                           'anchor', {'n','n'}, ...
                           'buffer',[0 40], ...
                           'nrow',1, ...
                           'fontsize',15);
    set(gca,'layer','top','XTickLabel',[])
    set(f,'CurrentAxes',ha(2))

    X = abs(Y)./repmat(sum(abs(Y),2),1,4) * 100;
    hold on
    bar(time_plot,...
        X,1,'stacked','EdgeColor','none','LineWidth',0.001);

    box on
    set(gca,'XMinorTick','on','YMinorTick','on',...
        'XTick',1:12:(size(Y,1)-1),...
        'XTickLabel',unique(time_dt.Year))
    ax = get(gca);
    ax.XAxis.MinorTickValues = 1:(size(X,1)-1+time_dt.Month(1));

    axis tight
    xlabel('Time')
    ylabel(sprintf('Contribution to \n density change (%%)'))
    set(gca,'layer','top') %,'LineWidth',2)

    % columnlegend(4,{'Sublimation','Dry compaction', 'Melt','Precipitation'},...
    %     'Location','SouthOutside');

    print(f, sprintf('%s/densification_driver_month_%s',c.OutputFolder,c.station),'-dtiff')
    if strcmp(vis,'off')
        close(f);
    end
end

%% plot year
if time_mod(end)-time_mod(1) >365
    time_start = datetime(datestr(time_mod(1)));  
    time_end = datetime(datestr(time_mod(end)));  
    time_new = datenum(time_start);
        i = 1;
        while time_new(end) < datenum(time_end)
           time_new(i) = datenum(time_start.Year + i-1, 6, 1);
           i = i + 1;
        end
    time_dt = datetime(datestr(time_mod));
    % ind_binned = time_dt.Month-time_dt.Month(1)+1 +12*(time_dt.Year - time_dt.Year(1));
    ind_binned = discretize(time_mod,time_new);
    ind_nan =isnan(ind_binned);

    Y = NaN(max(ind_binned),4);
    count =1;
    for i = [4 1 2 3]
        Y(:,count) = accumarray(ind_binned(~ind_nan), ...
            table2array(part_dens_hour(~ind_nan,i+1)));
        disp(part_dens_hour.Properties.VariableNames{i+1})
        disp(mean(Y(:,count)))

        count = count +1;
    end
    delta_rho_sum = accumarray(ind_binned(~ind_nan),part_dens_hour.check(~ind_nan));
    Z = [time_new(1:end-1)' Y];
    dlmwrite(sprintf('%s/yearly_densification_drivers.txt',c.OutputFolder),Z)

    if time_start>datetime(time_start.Year,6,10)
        Y(1,:) = NaN;
        delta_rho_sum(1) = NaN;
    end
    if time_end<datetime(time_end.Year,5,20)
        Y(end,:) = NaN;
        delta_rho_sum(end) = NaN;
    end

    Y1 = NaN(size(Y));
    Y1(Y>0) = Y(Y>0);
    Y2 = NaN(size(Y));
    Y2(Y<=0) = Y(Y<=0);
    Y_mean = repmat(nanmean(Y,1),size(Y,1),1);
    ind = Y >= Y_mean;
    ind(:,end) = ~ind(:,end);
    Ybis = Y;
    Ybis(ind) = 0;
    Y1bis = zeros(size(Ybis));
    Y1bis(Ybis>0) = Ybis(Ybis>0);
    Y2bis = zeros(size(Ybis));
    Y2bis(Ybis<=0) = Ybis(Ybis<=0);

    time_plot = datenum(time_dt.Year(1):(time_dt.Year(1)+size(Y,1)-1),6,1);

%     f = figure('Visible', vis);
%     ha = tight_subplot(2,1,0.03,0.2,0.1);
%     set(f,'CurrentAxes',ha(1))
    col = flipud(linspecer(4));
    hold on


    h4(1) = bar(time_plot+datenum(0.5,1,1),...
        Y1bis(:,1),1.001,'stacked','EdgeColor','none','LineWidth',0.001);
    for k = 1:3
        temp  = bar(time_plot+datenum(0.5,1,1),...
            [Y1(:,1:k), Y1bis(:,1+k)],1.001,'stacked','EdgeColor','none','LineWidth',0.001);
        h4(1+k) = temp(end);
    end

    h1 = bar(time_plot+datenum(0.5,1,1),...
        Y1,1.001,'stacked','EdgeColor','none','LineWidth',0.001);

 
    h3bis = bar(time_plot+datenum(0.5,1,1),...
        0*Y2,1.001,'stacked','EdgeColor','none','LineWidth',0.001);    
    h3ter = bar(time_plot+datenum(0.5,1,1),...
        0*Y2,1.001,'stacked','EdgeColor','none','LineWidth',0.001);
    h3 = bar(time_plot+datenum(0.5,1,1),...
        Y2,1.001,'stacked','EdgeColor','none','LineWidth',0.001);   
    uistack(h4,'top');

    h5 = bar(time_plot+datenum(0.5,1,1),...
        Y2bis,1.001,'stacked','EdgeColor','none','LineWidth',0.001);
    for i=1:4
        h1(i).FaceColor = col(i,:);
        h3(i).FaceColor = col(i,:);
        h3bis(i).FaceColor =[0.5 0.5 0.5];
        h4(i).FaceColor = h1(i).FaceColor +0.4*([1 1 1] - h1(i).FaceColor);
        h3ter(i).FaceColor = h3bis(i).FaceColor +0.4*([1 1 1] - h3bis(i).FaceColor);
        h5(i).FaceColor = h3(i).FaceColor +0.4*([1 1 1] - h3(i).FaceColor);
    end
    h_ghost = bar(time_plot+datenum(0.5,1,1),...
    Y1bis(:,1),1.001,'stacked','EdgeColor','none','LineWidth',0.001);
    h_ghost.FaceColor = [1 1 1];
    uistack(h_ghost,'bottom');

    h2 = stairs([time_plot time_plot(end)+365]-0.5,...
        delta_rho_sum([1:end end]),'k','LineWidth',2);

    plot([time_plot(1) time_plot(end)+365],[0 0], '--k')
    box on
    % 'XTick',1:12:(size(Y,1)-1),...
    %     'XTickLabel',unique(time_dt.Year))
    % ax = get(gca);
    %     ax.XAxis.MinorTickValues = 1:(size(Y,1)-1+time_dt.Month(1));
    set_monthly_tick(time_plot); 
    set(gca,'layer','top')
    set(gca,'XTickLabelRotation',45)

    axis tight
    xlim([time_plot(2) time_plot(end)+365])
    xlabel('Year')
    ylabel(sprintf('Contribution to \n density change (kg m^{-3})'),'Interpreter','tex')
    h_legend = legendflex([h1([4 3 2 1]), h2,h_ghost, h3bis(3), h3ter(3)],{'Precipitation', 'Melt-Refreeze','Grain-scale compaction',...
        'Sublimation','Total change',' ','above average','below average'}, ...
        'ref', gca,'box','off', ...
                           'anchor', {'n','n'}, ...
                           'buffer',[0 90], ...
                           'nrow',2, ...
                           'fontsize',15);

%     set(f,'CurrentAxes',ha(2))
%     set(gca,'Visible','off')
    % 
    % X = abs(Y)./repmat(sum(abs(Y),2),1,4) * 100;
    % hold on
    % bar(time_plot+datenum(0.5,1,1),...
    %     X,1.01,'stacked','EdgeColor','none','LineWidth',0.001);
    % 
    % box on
    % set_monthly_tick(time_plot); 
    % set(gca,'XTickLabelRotation',45)
    % axis tight
    % xlabel('Time')
    % ylabel(sprintf('Contribution to \n density change (%%)'))
    % set(gca,'layer','top') %,'LineWidth',2)

    % columnlegend(4,{'Sublimation','Dry compaction', 'Melt','Precipitation'},...
    %     'Location','SouthOutside');

%     print(f, sprintf('%s/densification_driver_year_%s',c.OutputFolder,c.station),'-dtiff')
%     if strcmp(vis,'off')
%         close(f);
%     end
end

%% plotting pore space

f = figure('Visible',vis);
[ha, l1,l2] = plotyy(time_mod,rho_avg,time_mod,pore_space);
hold(ha(1));
hold(ha(2));
if length(time_stamp)>1
scatter(ha(1), time_stamp, rho_stamp,120,'MarkerFaceColor','b','MarkerEdgeColor','b')
scatter(ha(2), time_stamp, ps_stamp,120,'MarkerFaceColor','r','MarkerEdgeColor','r')
end
set(l1,'LineWidth',2);
set(l2,'LineWidth',2);
ylabel(ha(1),'Average density (kg/m^3)')
ylabel(ha(2),'Pore space (m^3/m^2)')
xlabel('Date')
datetick('x','yyyy','keepticks')
box off
title('for the upper 20m')
set(ha(1),'XMinorTick','on','YMinorTick','on','XLim',time_mod([1 end]))
set(ha(2),'XMinorTick','on','YMinorTick','on','XLim',time_mod([1 end]))
print(f, sprintf('%s/pore_space_%s',c.OutputFolder,c.station), '-dtiff')
    close(f);



end
