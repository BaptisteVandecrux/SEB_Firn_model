function [rho_avg] = Calculate20mAvgDensity(snowc, snic, slwc, snowbkt, rho, c)
    thickness_weq = snowc + snic +slwc;
    thickness_weq = [snowbkt'; thickness_weq];

    thickness_act = snowc.*(c.rho_water./rho) + snic.*(c.rho_water./c.rho_ice);
    thickness_act =[snowbkt'*c.rho_water/315; thickness_act];

    depth_act = cumsum(thickness_act );
    depth_weq = cumsum(thickness_weq);

    rho_all= (snowc + snic)./...
                (snowc./rho + snic./c.rho_ice);
    rho_all = [315*ones(size(snowbkt')); rho_all];

    depth_alt = [depth_act; 20];
    rho_alt = [rho_all; NaN];

    [depth_alt, ind_sorted] = sort(depth_alt,1);
    rho_alt = rho_alt(ind_sorted);

    % we then interpolate at each time step
    ind_20m = isnan(rho_alt);
    ind_next =  boolean([0; ind_20m(1:end-1,:)]);

try    rho_alt(ind_20m) = rho_alt(ind_next);
catch me
    khdfhv= 0;
end

    thickness_alt = depth_alt;
    thickness_alt(2:end,:) = depth_alt(2:end,:) - depth_alt(1:end-1,:);

    rho_temp = NaN(size(rho_alt));
    rho_temp(depth_alt <= 20) = rho_alt(depth_alt(:) <= 20);
    rho_avg = nansum(rho_temp .* thickness_alt,1) ./20;
end