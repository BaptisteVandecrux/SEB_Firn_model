function [Tsurf_reset, T_ice_reset] = ...
    ResetTemp(depth_thermistor, LRin, LRout, T_obs, rho, T_ice, time, k, c)

% ========== resets surface temperature ================
if ~isnan(LRout(k)) && ~isnan(LRin(k))
    Tsurf_reset = ((LRout(k) - (1-c.em)*LRin(k)) /(c.em*c.sigma))^(1/4);
else
    Tsurf_reset = NaN;
end

% ================= resets temperature profile =========================

    % if there is thermistor record for the first time step, then reads 
    % initial subsurface conditions from AWS data
    T_ice_reset = NaN(c.jpgrnd,1);
    if sum(~isnan(T_obs(k,:)))>1
        depth = depth_thermistor(k,(depth_thermistor(k,:)~=0))';
        oldtemp = T_obs(k,(depth_thermistor(k,:)~=0))';
    end
    
    % calculates the new depth scale in mweq
    depth_weq = cumsum(c.cdel);

    % calculates the new depth scale in real m
%     depth_act = depth_weq .*c.rho_water ./rho(:,k);
    depth_act = cumsum(c.cdel .*c.rho_water ./rho(:,k));

    % Here we add an initial surface temperature
    depth_act = [0; depth_act];
    depth_weq = [0; depth_weq];
    if depth(1) ~= 0
        %if the input density profile does not contain surface temperature,
        %then we use Tsurf_in to initiate the temperature profile
        depth = [0; depth];
        oldtemp = [Tsurf_reset - c.T_0; oldtemp];
    end
    % the old scale is converted from m to mweq by interpolation
    oldscale_weq = interp1(depth_act,depth_weq,depth);

    % we limit the new scale to the values that are given within the oldscale
    newscale_weq = depth_weq(depth_weq<=oldscale_weq(end));

    % the temperature is interpolated on each stamp of the new scale
    newtemp = interp1(oldscale_weq,oldtemp,newscale_weq);
    newtemp = newtemp + c.T_0; %going back to K

    % giving the observed temperature (on appropriate scale) as initial value
    % for the subsurface temperature profile
    % There might be a shift to apply depending on whether the first value in
    % subsurface column represent the temp at depth 0 (=Tsurf) or at depth 
    % c.cdel(1). To be asked.
    T_ice_reset(1:length(newtemp)) = newtemp;

    % fills the rest of the temperature profile withan interpolation that
    % ends with Tdeep at the bottom of the profile

    if length(newtemp)<length(c.cdel)
        d1 = length(newtemp);
        T1 = newtemp(d1);
        d2 = c.jpgrnd;
        T2 = c.Tdeep_AWS + c.T_0;
        ind = length(newtemp)+1:(c.jpgrnd-1);
        T_ice_reset(length(newtemp)+1:(c.jpgrnd-1)) = ...
            (T1-T2)/(d2-d1)^2*(d2-ind).^2 + T2;
        T_ice_reset(c.jpgrnd) = c.Tdeep_AWS+ c.T_0;

%         T_ice_reset(length(newtemp)+1:(c.jpgrnd-1)) = NaN;
    end
    
%     figure
%     scatter(oldscale_weq,oldtemp,'o')
%     hold on
%     stairs(depth_weq(1:end-1),T_ice(:,k,1)-c.T_0, 'LineWidth',2)
%     stairs(depth_weq(1:end-1),T_ice_reset-c.T_0, 'LineWidth',2)
%     legend('data','before reset','after reset','Location','South')
%     xlabel('Depth (m weq)')
%     ylabel('Temperature (deg C)')
%     title(sprintf('2/ %s',datestr(datenum(time(k),0,0))))
%     view([90 90])
    
% removing non-freezing temperatures (just in case)
subsurfmelt = find(T_ice_reset(:) > c.T_0);
if sum(subsurfmelt )>0
    T_ice_reset(subsurfmelt) = c.T_0;
end
%% Comes from outside the function just after it
%                     figure
                    depth_act = cumsum(c.cdel .*c.rho_water ./rho(:,k));
                    depth_act = [0; depth_act];

%                     scatter(depth_thermistor(k,depth_thermistor(k,:)~=0),...
%                         T_ice_obs(k,depth_thermistor(k,:) ~= 0), 'o')
%                     hold on
%                     stairs(depth_act(1:end-1),T_ice(:,k,j)-c.T_0)
                    
                T_ice(~isnan(T_reset),k,j) = T_reset(~isnan(T_reset));
%                     stairs(depth_act(1:end-1),T_ice(:,k,j)-c.T_0)
%                     legend('data','before reset','after reset','Location','South')
%                     xlabel('Depth (m)')
%                     ylabel('Temperature (deg C)')
%                     title(sprintf('%s',datestr(datenum(time(k),0,0))))
%                     view([90 90])
   
                    [zso_capa, zso_cond] = ice_heats (c);
                    [grndc, grndd, ~, ~]...
                        = update_tempdiff_params (rho(:,k), Tdeep(j)                    ...
                        , snowc, snic, T_ice(:,k,j), zso_cond, zso_capa, c);

end
