function [melt_mweq, sublimation_mweq, snowthick] = ...
    MassBudget (meltflux, H_snow, H_rain, LHF, ...
    snowfall, rainfall, snowthick, elev, j, k, c)
% UNUSED

% MassBudget: Calculates the surface mass budget including diverse surface 
% sources and sinks of liquid water, snow and ice.
%
% Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================
% in mweq
    melt_mweq = -meltflux(k,j)*c.dt_obs/c.dev/c.L_fus/c.rho_water;   
    sublimation_mweq = LHF(k,j)*c.dt_obs/c.dev/c.L_sub/c.rho_water; % in mweq
    % positive LHF -> deposition -> dH_subl positive


end
