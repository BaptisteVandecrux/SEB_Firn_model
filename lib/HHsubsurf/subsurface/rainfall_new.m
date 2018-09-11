function [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil ] ...
    = rainfall_new (zraind, psnowc, psnic, pslwc, pdgrain ...
    , prhofirn, ptsoil, pts, c)

%add_rainfall: Routine that adds rain water input on top of the 
%subsurface column. Update BV2017: No mass shift anymore. Water is just
%added to the water fraction of the top layer and the total mass of the layer is
%allowed to change.
%	
%   input:
%         zraind - amount (m weq) of input rain.
%         
%         prhofirn - vector of firn (snow) density (kg/m^3) 
%
%         psnowc, psnic, pslwc - vectors of respectively snow, ice and
%         liquid water part for each subsurface layer (m weq).
%
%         ptsoil - vector of subsurface temperature (K)
%
%         pdgrain - Vector of layer grain size. see graingrowth function.
%
%         pts - surface temperature (K)
%
%         c - Structure containing all the physical, site-dependant or user
%         defined variables.
%
%   output:
%          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil]
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================

if zraind > c.smallno
    T_rain_in = max(pts,c.T_0);
    ptsoil(1) = ((psnowc(1) + psnic(1) + pslwc(1)) * ptsoil(1)  + zraind * T_rain_in) ...
        / (psnowc(1) + psnic(1) + pslwc(1) + zraind);

    pslwc(1) = pslwc(1) + zraind;
    % rain switch the snowbkt into the first layer to mimic the fresh snow
    % getting wet

end