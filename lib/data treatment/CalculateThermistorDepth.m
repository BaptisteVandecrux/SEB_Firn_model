function [depth_therm_modscale, T_ice_obs] = CalculateThermistorDepth(time_mod,...
    depth_mod, H_surf,compaction, depth_obs, T_ice_obs, station_name)
% This function calculates the thermistor depth based on a given surface
% height record (and therefore the burial rate), a maintenance report 
% and a compaction grid which indicates the compaction occuring between the
% buried sensores. Eventually removes the observed subsurface temperature
% that are less than 0.5 m below the surface.
%   
% Inputs:
% - time_mode: vector containing modelled time stamps
% - depth_absolute: matrix containing, for each model layer and at each
% time step, the absolute depth of the layer. The depth scale is anchored 
% at the bottom of the model, at a depth below which no firn-related
% compaction should occur. The zero of this scale is located at the initial
% surface level at the beginning of the simulation. Depth increase
% downward.
% - H_surf: Modelled surface height
% - depth_obs: Previous (rougher) estimation of sensors' depth (should be
% removed at some point).
% - station_name: Used to locate the maintenance file.
% - T_ice_obs: Observed subsurface temperature
%     
% Ouputs:
% - depth_therm_modscale: matrix containing the depth of the thermistors on
% the model's absolute depth scale (see depth_absolute description)
%
% 16-10-2018
% B. Vandecrux
% b.vandecrux@gmail.com
% ========================================================================
    
step = 12; %for faster plotting

    if strcmp(station_name,'CEN_THM')
        depth_therm_modscale = NaN(50,length(time_mod));
        maintenance = table();
        maintenance.date = datenum('26-Jul-2017 15:00:00');
        for i =1:50
            maintenance.(['NewDepth', num2str(i), 'm'])= ...
                depth_obs(i,find(~isnan(depth_obs(i,:)),1,'first'));
        end
    else
        maintenance = ImportMaintenanceFile(station_name);
        maintenance=standardizeMissing(maintenance,999);
    end
    ind_therm_fields = find(contains(maintenance.Properties.VariableNames,'Depth'));
    ind_no_new_depth = isnan(nanmean(table2array(...
        maintenance(:,ind_therm_fields)),2));
    maintenance(ind_no_new_depth,:) = [];

    varnames = maintenance.Properties.VariableNames;
    depth_therm_modscale = NaN(size(depth_obs,1),length(time_mod)); 
    for j = 1:size(maintenance,1)
        ind_time = find(time_mod >= maintenance.date(j));                
        [~, ind_start] = min(abs(time_mod - maintenance.date(j)));
        if j<size(maintenance,1)
            if maintenance.date(j+1)<time_mod(1)
                continue
            end
            [~, ind_end] = min(abs(time_mod - maintenance.date(j+1)));
            ind_end=ind_end-1;
        else
            ind_end=length(time_mod);
        end
        disp(datestr(time_mod([ind_start ind_end])))

        % assigning the depth of the thermistors on the model's absolute
        % depth scale
        if time_mod(1) > maintenance.date(j)+1           
            ind = find(~isnan(depth_obs(1,:)),1,'first');
            depth_therm_modscale(:,ind_start) = ...
                depth_obs(:,ind);
        else
            depth_therm_modscale(:,ind_start) = ...
                table2array(maintenance(j,ind_therm_fields));
        end
        % for each time step following the installation
        for i = (ind_start+step):step:ind_end
            % we first assume that the sensor is still buried at the
            % same absolute depth. Since we work with absolute depth
            % the material accumulating or sublimating at the surface
            % do not interfer
            depth_therm_modscale(:,i) = max(0,...
                depth_therm_modscale(:,i-step) + (H_surf(i)-H_surf(i-step)));
   % however the material between the sensor and the bottom of
            % the model was still subject to compaction within that
            % time step which should decrease the absolute depth of the
            % sensor
%             if max(depth_mod(:,i))<max(depth_therm_modscale(:,i))
%                 depth_mod_tmp=[depth_mod(:,i); max(depth_therm_modscale(:,i))];
%                 compaction_tmp=[compaction(:,i); compaction(end,i)];
%             else
                depth_mod_tmp = depth_mod(:,i);
                compaction_tmp = compaction(:,i-step);
%             end
                                        
            temp = [depth_therm_modscale(:,i); depth_mod_tmp];
            [new_depth,ind_depth] = sort(temp);
            for d = 1:length(depth_therm_modscale(:,i))
                ind_therm_in_new_depth(d) = find(ind_depth == d);
            end
            % here "ConvertToGivenDepthScale" works with depth interval
            % starting at 0 positive downward, so we feed that function
            % with the model depth and desired depth relative to the
            % surface
            compaction_2 = ConvertToGivenDepthScale(depth_mod_tmp,...
                compaction_tmp, new_depth, 'extensive')*step;

            % now we update our first guess by adding the sum of all
            % compaction that occured bellow the sensor
            for kk = 1:size(depth_obs,1)
                depth_therm_modscale(kk, i) = depth_therm_modscale(kk, i)...
                    + sum(compaction_2(ind_therm_in_new_depth(kk)+1:end));
            end
        end

    end
    for kk = 1:size(depth_obs,1)            
        % interpolating between the steps
        ind_nan = find(isnan(depth_therm_modscale(kk, :)));
        ind_no_nan = find(~isnan(depth_therm_modscale(kk, :)));
         depth_therm_modscale(kk,ind_nan)=interp1(ind_no_nan,...
            depth_therm_modscale(kk, ind_no_nan),...
            ind_nan,'linear');
    end
    % sorting the depth 
    for i = 1:length(time_mod)
        [depth_therm_modscale(:,i), ind]= ...
            sort(depth_therm_modscale(:,i),'ascend');
        T_ice_obs(:,i) = T_ice_obs(ind,i);
    end
end
    