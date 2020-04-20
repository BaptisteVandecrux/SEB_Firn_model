function  [SRout_mdl, SRnet, T_ice, meltflux_internal, dH_melt_internal] = ...
    SRbalance (SRout_mdl, SRin, SRnet,...
    z_icehorizon,snowthick, T_ice, rho, j, k, c)

% SRbalance: Calculates the amount of Shortwave Radiation that is 
% penetrating at each layer (SRnet).
% uses it to warm each layer and eventually calculates the melt that is
% produced by this warming
%
% Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================
%extinction coefficient of ice 0.6 to 1.5 m-1
%extinction coefficient of snow 4 to 40 m-1
% Greufell and Maykut, 1977, Journal of Glaciology, Vol. 18, No. 80, 1977

% %radiation absorption in snow
% SRnet(snow_layer) = (SRin - SRout_mdl)...
%     *exp(-ext_snow*depth(snow_layer));
% 
% %radiation absorption in ice layers underneath the snowpack
%  SRnet(ice_layer) = (SRin-SRout).*...
%         exp(-ext_snow*snowthick).*...
%         exp(-ext_ice*(depth(ice_layer) - snowthick));

% NOT USED AT THE MOMENT
% should be updated for non regular layer thickness
    
if c.elev_bins ~= 1 && k > 1
    % in a transect, SRout is calculated using SRin (calculated from SRin_AWS)
    % and fixed values for albedo for snow and ice
    if snowthick(k,j) > 0
        SRout_mdl(k,j) = c.alb_snow*SRin(k,j) ;
    else
        SRout_mdl(k,j) = c.alb_ice*SRin(k,j);
    end
end
%radiation absorption in snow

SRnet(1:(z_icehorizon+1)) = (SRin(k,j) - SRout_mdl(k,j))...
    *exp(-c.ext_snow*(0:z_icehorizon)*c.dz_ice);

%radiation absorption in underlying ice
if z_icehorizon < c.z_ice_max
    SRnet((z_icehorizon+2):c.z_ice_max+1) = (SRin(k,j)-SRout_mdl(k,j)).*...
        exp(-c.ext_snow*snowthick(k,j)).*...
        exp(-c.ext_ice*(((z_icehorizon+2):(c.z_ice_max+1))*...
        c.dz_ice - snowthick(k,j)));
end
% specific heat of ice 
%(perhaps a slight overestimation for near-melt ice temperatures (max 48 J/kg/K))
if k==1
    c_i = 152.456 + 7.122 * T_ice(:,1,j) ;   
%snow & ice temperature rise due to shortwave radiation absorption
else
    % Specific heat of ice (a slight overestimation for near-melt T (max 48 J kg-1 K-1))
    c_i = 152.456 + 7.122 * T_ice(:,k-1,j);
    
    T_ice(1:c.z_ice_max,k,j) = T_ice(1:c.z_ice_max,k-1,j) +...
        c.dt_obs./c.dev./rho(1:c.z_ice_max,k)./c_i(1:c.z_ice_max).*...
        (SRnet(1:c.z_ice_max)-SRnet(2:(c.z_ice_max+1)))./c.dz_ice;
    
%      a=...
%         (SRnet(1:c.z_ice_max)-SRnet(2:(c.z_ice_max+1)))./c.dz_ice;
    
    T_ice(c.z_ice_max+1,k,j) = T_ice(c.z_ice_max+1,1,j);

end

% finding where/how much melt occurs
subsurfmelt = (T_ice(:,k,j) > c.T_0);
nosubsurfmelt = (T_ice(:,k,j) <= c.T_0);
dT_ice = T_ice(:,k,j) - c.T_0;
dT_ice(nosubsurfmelt) = 0;

meltflux_internal_temp = rho(:,k) .* c_i .* dT_ice / c.dt_obs * c.dev * c.dz_ice;
meltflux_internal = sum(meltflux_internal_temp(1:c.z_ice_max));

dH_melt_internal = zeros (c.jpgrnd,1);
dH_melt_internal(:) = -c_i.*dT_ice./c.L_fus * c.dz_ice;
dH_melt_internal(1) = 0;

if sum(subsurfmelt) > 0
    T_ice(subsurfmelt,k,j) = c.T_0;   % removing non-freezing temperatures
end
% Reduce sub-surface density due to melt? Will most likely cause model instability...      
end
