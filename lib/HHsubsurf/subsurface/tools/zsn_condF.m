% PLA densification (Feb 2015
function [zsn_condF] = zsn_condF(prhofirn, c)
% The F in the name zsn_condF stands for "function"

% update BV2017
zsn_condF = 2.22362*(prhofirn./1000).^1.88;
% zsn_condF2 = cpice(0)*c.rho_ice*c.zdifiz*((prhofirn/c.rho_water).^1.88);
zsn_condF2 = (1-(c.rho_ice - prhofirn)/c.rho_ice)*2.1;

% snow thermal conductivity [J/s/m/K]
% Yen, Y. (1981), Review of thermal properties of snow, ice and sea ice,
% Rep. 81-10, U. S. Army Cold Reg. Res. and Eng. Lab. (CRREL), Hanover, N. H
end
