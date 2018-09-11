function [snowthick, z_icehorizon] = ...
    UpdateSnowThickness(snowthick,z_icehorizon, k, j, c)

% UpdateSnowThickness: Calculates variation of snow thickness since initial 
% conditions (snowthick_AWS).
%
% Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================
  
if k==1
    % Initial snowpack thickness
    %Initial snow depth
    if c.elev_bins ==  1
        snowthick(1,j) = c.snowthick_ini ;
    elseif elev(j) < c.ELA
        snowthick(1,j) = c.snowthick_ini + c.gradsnowthick*(elev(j)-c.ELA);
    end
    %Snow to ice transition
    z_icehorizon = min(c.z_ice_max, floor(snowthick(1,j)/c.dz_ice));
end
% Update BV2017: sensor height from input data
if k>1      
    snowthick(k,j) = snowthick(k-1,j);       % will be updated at end of time loop
end
    
end
