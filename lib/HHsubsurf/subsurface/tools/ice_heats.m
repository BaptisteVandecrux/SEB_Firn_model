function [zso_capa, zso_cond] = ice_heats (ptsoil, c)
%ice_heats:computes the subsurface thermal capacity [J/K] and thermal 
% conductivity zso_cond [J/S/M/K] from the subsurface temperature diffusivity
% c.zdifiz [M**2/S]
% Input:
%   ptsoil - subsufrace temperature in Kelvin
%   c - structure containing the constants
% Output:
%   zso_capa - Volumetric heat capacity of ice for each of the layers in
%   J/kg/K. Note that it is the same for snow!
%   zso_cond - Thermal conductivity of ice for each of the layer, in W/m/K
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================

% Update BV2017
% Volumetric heat capacity of ice. In J/m^3/K.
zso_capa = zsn_capaF(c.rho_ice, ptsoil);

% update BV2017
% thermal conductivity of pure ice in W/m/K according to Yen 1981
zso_cond = 9.828*exp(-0.0057*ptsoil);

% old version
% zso_cond = zso_capa*c.zdifiz;

end
