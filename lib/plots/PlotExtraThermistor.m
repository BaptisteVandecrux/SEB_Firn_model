function [] = PlotExtraThermistor(path_list,time_mod,depth_obs_s,...
    depth_obs_m,T_ice_obs, T_ice_mod, H_surf,station, OutputFolder,...
    T_subsurf_mod, depth_act,vis)


%% Ploting observed temperature
 vis = 'on';
 f=figure('Visible',vis);%('outerposition',[1 -1  25 25]);
 set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,[0.02 0.02],[0.25 0.15],[0.07 0.1]);
     load(strcat(path_list{1},'/run_param.mat'))

     scale = 'model';
%      scale = 'station';
ylimit = 10;
    step = 72;
    ind_sta = [1 2 8];
 for count = 1:length(ha)
     tmp =mod(count,3);
     if tmp == 0
         tmp = 3;
     end
     ii = ind_sta(tmp);

    switch ii
        case 1 
            xlimit = datenum(2007,[6 12],1);
            ylimit_plot=[-8 3];
            title_text='a) CP1';
        case 2 
            xlimit = datenum(2016,[4 12],1);
                        ylimit_plot=[-8 3];
            title_text='b) DYE-2';
        case 3
            xlimit = datenum(2013,[4 12],1);
            ylimit_plot=[-8 3];
            title_text='c) Summit';
    end
    T_ice_obs{ii}(:,or(time_mod{ii}<xlimit(1),time_mod{ii}>xlimit(2))) = NaN;

    set(f,'CurrentAxes',ha(count))



switch count
    case {1, 2, 3}
        title(title_text)
        TT = ones(size(depth_obs_s{ii},1),1) * time_mod{ii}';
        [~, ind_bot] = min(abs(...
            depth_obs_s{ii}(:,~or(time_mod{ii}<xlimit(1),time_mod{ii}>xlimit(2)))...
            - (ylimit+1)));
        ind_lim = min(max(ind_bot)+1, size( T_ice_obs{ii},1));

        switch scale
            case 'model'
                depth_temp = depth_obs_m{ii};

            case 'station'
                depth_temp = depth_obs_s{ii};
                for i = 1: size(depth_temp,1)
                    depth_temp(i,:) = depth_temp(i,:) - Surface_Height{ii}';
                end
        end

        temp_plot =    T_ice_obs{ii}(1:ind_lim, 1:step:end);

    case {4, 5 ,6}

        TT = ones(size(depth_act{ii},1),1) * time_mod{ii}';
        [~, ind_bot] = min(abs(depth_act{ii} - (ylimit+1)));
        ind_lim = min(max(ind_bot)+1, size( T_ice_mod{ii},1));
        step = 72;
        depth_temp = depth_act{ii};
       temp_plot = T_ice_mod{ii}(1:ind_lim, 1:step:end)-273.15;
        
    case {7, 8, 9}
        TT = ones(size(depth_obs_s{ii},1),1) * time_mod{ii}';
        [~, ind_bot] = min(abs(...
            depth_obs_s{ii}(:,~or(time_mod{ii}<xlimit(1),time_mod{ii}>xlimit(2)))...
            - (ylimit+1)));
        ind_lim = min(max(ind_bot)+1, size( T_ice_obs{ii},1));

        switch scale
            case 'model'
                depth_temp = depth_obs_m{ii};

            case 'station'
                depth_temp = depth_obs_s{ii};
                for i = 1: size(depth_temp,1)
                    depth_temp(i,:) = depth_temp(i,:) - Surface_Height{ii}';
                end
        end

        T_diff = T_subsurf_mod{ii} - T_ice_obs{ii};
        temp_plot =   T_diff(1:ind_lim, 1:step:end);        
end
    col = PlotTemp(TT(1:ind_lim, 1:step:end),...
        depth_temp(1:ind_lim, 1:step:end),...
       temp_plot,...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'ShowLegend','no',...
        'cmap','jet',...
        'Interp','on',...
        'XLabel','',...
        'YLabel',' ',...
        'CLabel','Firn temperature (^oC)',...
        'Range', -40:1:0,...
        'FlatSurface','no');
%     col.FontSize = 12;
    switch count
        case 1
        col.Position(1) = 0.92;
        col.Position(4) = 0.3 ;
        col.Position(2) = 0.6;

        case 7
            for i = 1:6
            freezeColors(ha(i))
            end
            colormap(ha(count),'parula')
        col.Position(1) = 0.92;
        col.Position(4) = 0.3 ;
        col.Position(2) = 0.08;
        otherwise
        col.Position(1) = col.Position(1)+ 3;
    end
    
    switch scale
        case 'model'
            plot(time_mod{ii}, - H_surf{ii},'LineWidth',2)
            plot(time_mod{ii}, 10 - H_surf{ii},':','LineWidth',2)
        case 'station'
            plot(time_mod{ii}, - Surface_Height{ii}','LineWidth',2)
            plot(time_mod{ii}, 10 - Surface_Height{ii}',':','LineWidth',2)
    end
    
%     col.TickLength=col.TickLength*5;
    temp = col.YTickLabel;
    for k = 1:length(temp)
        if (k-1)/5==floor((k-1)/5)
        else
            temp(k,:)=' ';
        end
    end
    col.YTickLabel=temp;

    ylim(ylimit_plot)
    xlim(xlimit)


    aux = datevec(time_mod{ii});
    ind_in = ~or(time_mod{ii}<xlimit(1),time_mod{ii}>xlimit(2));
    aux = aux(ind_in,:);
    aux(:,3) = 1;
    aux(:,4:end) = 0;
    aux = unique(aux,'rows');
    set(gca,'XTick',datenum(aux))

    datetick('x','mm-yyyy','keeplimits','keepticks')
    if count<7
        set(gca,'XTickLabel','')
    end
end

print(f, sprintf('%s/T_ice_obs_extra',OutputFolder), '-dpng')
    if strcmp(vis,'off')
        close(f)
    end
end
    