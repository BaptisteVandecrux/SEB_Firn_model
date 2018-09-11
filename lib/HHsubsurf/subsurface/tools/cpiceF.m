function [cpice] = cpiceF(ptsoil)
%cpice: right now constant heat capacity of ice. Can be modified to be
%dependant on ice caracteristics.
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================
% Specific heat capacity of ice J/kg/K
% cpice = 2.108e+03  ;    % standard value at 0degC 

% Update BV2017
cpice = (2.7442 + 0.1282 .* ptsoil)/18*1000;
% Yen 1981 originally in J/mol/K converted to J/kg/K
% 1kg of ice has same capacity as 1kg snow.

end
