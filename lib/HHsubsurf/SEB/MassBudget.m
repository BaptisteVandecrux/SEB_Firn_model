function [dH_melt_weq, H_melt_weq, dH_subl_weq, H_subl, H_snow, H_rain, snowthick] = ...
    MassBudget (meltflux, H_subl, H_melt_weq, H_snow, H_rain, LHF, ...
    snowfall, rainfall, snowthick, elev, j, k, c)
% MassBudget: Calculates the surface mass budget including diverse surface 
% sources and sinks of liquid water, snow and ice.
%
% Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================
% in mweq
dH_melt_weq = -meltflux(k,j)*c.dt_obs/c.dev/c.L_fus/c.rho_water;   

if k > 1
    H_melt_weq(k,j) = H_melt_weq(k-1,j) + dH_melt_weq;  %in mweq
end

dH_subl_weq = LHF(k,j)*c.dt_obs/c.dev/c.L_sub/c.rho_water; % in mweq
% positive LHF -> deposition -> dH_subl positive
if k > 1
    H_subl(k,j) = H_subl(k-1,j) + dH_subl_weq;
end

% Total precipitation in m of water
if elev(j) <c.prec_cutoff
    snowfall(k,j) = snowfall(k,j)/2;
    rainfall(k,j) = rainfall(k,j)/2;
end

if k > 1
    %H_snow in m of snow
    H_snow(k,j) = H_snow(k-1,j) + snowfall(k,j)/c.rho_snow(k,j)*c.rho_water;
else
    H_snow(k,j) = 0;
end
if k > 1
    %in mweq
    H_rain(k,j) = H_rain(k-1,j) + rainfall(k,j);
end

% UPDATE BV 2016: Now the surface height change is calculated after
% subsurface scheme
   
if k > 1
    %in mweq
    snowthick(k,j) = snowthick(k-1,j) + ...
        max ( 0, snowfall(k,j)*c.rho_snow(k,j)/c.rho_water + dH_melt_weq + dH_subl_weq);
end
end
