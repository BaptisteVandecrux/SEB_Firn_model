function [rho_avg, CC20] = Calculate20mAvgDensity(snowc, snic, slwc, snowbkt, rho,ptsoil, c)
    % function that calculates 20 m average density and total cold content 
    % to 20 m depth
    % Baptiste Vandecrux 04-08-2019
    %

    %% calculate depth scale and bulk density
    thickness_weq = snowc + snic +slwc;
    thickness_weq = [snowbkt'; thickness_weq];

    thickness_act = snowc.*(c.rho_water./rho) + snic.*(c.rho_water./c.rho_ice);
    thickness_act =[snowbkt'*c.rho_water/315; thickness_act];

    depth_act = cumsum(thickness_act );
    depth_weq = cumsum(thickness_weq);

    rho_all= (snowc + snic)./...
                (snowc./rho + snic./c.rho_ice);
    rho_all = [315*ones(size(snowbkt')); rho_all];

    %% first calculation as in Vandecrux et al., 2018
    depth_alt = [depth_act; 20];
    rho_alt = [rho_all; NaN];

    [depth_alt, ind_sorted] = sort(depth_alt,1);
    rho_alt = rho_alt(ind_sorted);

    % we then interpolate at each time step
    ind_20m = isnan(rho_alt);
    ind_next =  boolean([0; ind_20m(1:end-1,:)]);

    rho_alt(ind_20m) = rho_alt(ind_next);

    thickness_alt = depth_alt;
    thickness_alt(2:end,:) = depth_alt(2:end,:) - depth_alt(1:end-1,:);

    rho_temp = NaN(size(rho_alt));
    rho_temp(depth_alt <= 20) = rho_alt(depth_alt(:) <= 20);
    rho_avg = nansum(rho_temp .* thickness_alt,1) ./20;
    
    %% second calculation as in Vandecrux et al. 2019
    ColdContent = @(rho,thickness,T_firn) ...
    2.108 .* rho .* thickness .* (c.T_0 - T_firn);
%    kJ/kg/K * kg/m3 * m            * K
% in kJ/m2

    [depth_alt, rho_alt, rho_20] = ...
        InterpGridColumnWise(depth_act, rho_all, 20, 'next','VolumeWeighted',[]);
    [~, T_ice_alt, T_20] = ...
        InterpGridColumnWise(depth_act, [ptsoil(1); ptsoil], 20, 'lin','MassWeighted', rho_alt);

    CC20  = ColdContent(rho_20, 20*ones(size(rho_20)), T_20);
    CC20_bis  = sum(ColdContent(rho_alt(depth_alt<=20), ...
        thickness_alt(depth_alt<=20),...
        T_ice_alt(depth_alt<=20)));
    
    if rho_avg~=rho_20 || abs(CC20/1000 - CC20_bis/1000)>0.1
        disp([CC20/1000 CC20_bis/1000])
        error('Differing results')
    end
    

end