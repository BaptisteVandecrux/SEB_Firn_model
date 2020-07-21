function [] = SurfaceHeightBudget(path_list,time_mod, ...
    depth_act, depth_act_2, H_surf, station, ...
   compaction, Surface_Height,write_surface_height, OutputFolder, vis)
%% depth station
ylimit = 10;
 f=figure('Visible',vis);%('outerposition',[1 -1  25 25]);
 set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,0.02,[0.15 0.02],[0.07 0.1]);
for ii = 1:length(path_list)
            load(strcat(path_list{ii},'/run_param.mat'))
 
   depth_from_surf = depth_act_2{ii};
    [~, ind_bot] = min(abs( 4 - depth_from_surf(:,1)));

    depth_bot_station = [];
    depth_bot_station(1) = depth_act{ii}(ind_bot,1);
    for i = 1:length(time_mod{ii})-1
        depth_bot_station(i+1) = depth_bot_station(i) + sum(compaction{ii}(ind_bot:end, i));
        [~, ind_bot] = min(abs( depth_bot_station(i+1) - depth_act{ii}(:,i+1)));
    end

    % plotting
    set(f,'CurrentAxes',ha(ii))
    hold on
    plot(time_mod{ii}, H_surf{ii})
    plot(time_mod{ii},-depth_bot_station)

    set_plot_style(ii, station, time_mod, 'Height above initial surface level (m)');
end
	    print(f, sprintf('%s/depth_station',OutputFolder), '-dpng')
    if strcmp(vis,'off')
        close(f)
    end
   
    %% Surface Height budget
f=figure('Visible',vis);%('outerposition',[1 -1  25 25]);
 set(gcf,'Position',[0    1.0583   36.1421   17.2508])
 ha = tight_subplot(3,3,0.03,[0.15 0.15],[0.07 0.1]);
for ii = 1:length(path_list)    
    depth_from_surf = depth_act_2{ii};
    [~, ind_bot] = min(abs( 4 - depth_from_surf(:,1)));
    depth_bot_station = [];
    depth_bot_station(1) = depth_act{ii}(ind_bot,1);
    
    for i = 1:length(time_mod{ii})-1
        depth_bot_station(i+1) = depth_bot_station(i) + sum(compaction{ii}(ind_bot:end, i));
        [~, ind_bot] = min(abs( depth_bot_station(i+1) - depth_act{ii}(:,i+1)));
    end
      H_surf_comp = H_surf{ii}-H_surf{ii}(1) + depth_bot_station'-depth_bot_station(1);

    Surf_Height_cor = ...
        Surface_Height{ii}-Surface_Height{ii}(find(~isnan(Surface_Height{ii}),1,'first'));
    Surf_Height_cor = Surf_Height_cor - depth_bot_station'+depth_bot_station(1);
       
    if write_surface_height == 1
        data = {Surface_Height{ii}, Surf_Height_cor, depth_bot_station'};
        namefile = [OutputFolder '\Surface_Height_' station{ii} '.nc'];
        varname = {'Surface_Height_obs', 'Surface_Height_cor', 'station_depth_mod'};
        unit = {'m' 'm' 'm'};
        long_varname = {'Surface height observed by the station, positive upward, with reference to the bottom of the mast.',...
            'Surface height, positive upward, corrected so that the referential is at a depth where firn compaction does not occur anymore.',...
            'Modelled depth of the bottom of the station, positive downward.'};

        WriteNC_1D(namefile, time_mod{ii}, data, varname, unit, long_varname);
    end
    % plotting
    set(f,'CurrentAxes',ha(ii))
    hold on
    plot(time_mod{ii}, Surface_Height{ii}-Surface_Height{ii}(find(~isnan(Surface_Height{ii}),1,'first')),'Color','b','LineWidth',2);
    plot(time_mod{ii},H_surf{ii},'Color', 'k', 'LineWidth', 2);
    plot(time_mod{ii}, H_surf_comp,'Color', 'r', 'LineWidth', 2);

    axis tight
    ylim([min(min(H_surf{ii},Surface_Height{ii}-Surface_Height{ii}(find(~isnan(Surface_Height{ii}),1,'first'))))*1.1 ...
        max(max([H_surf{ii}, ...
        Surface_Height{ii} - Surface_Height{ii}(find(~isnan(Surface_Height{ii}),1,'first')),...
        H_surf_comp],[], 2))*1.1])

    if ii ==1
    legendflex({'Observed surface height',...
    'Modelled surface height (origin at bottom of column)', ...
    'Modelled surface height (origin at bottom of station mast)'},...
    'ref', gcf, ...
    'anchor', {'n','n'}, ...
    'buffer',[0 0], ...
    'ncol',2, ...
    'fontsize',15);
    end
    
    set_plot_style(ii, station, time_mod, 'Height above initial surface level (m)');

end
	    print(f, sprintf('%s/surface height budget',OutputFolder), '-dpng')
    if strcmp(vis,'off')
        close(f)
    end
    end