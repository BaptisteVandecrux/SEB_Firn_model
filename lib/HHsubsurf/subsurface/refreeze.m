function [psnic, pslwc, ptsoil, zrfrz]...
    = refreeze (psnowc, psnic, pslwc, ptsoil, c )
%refreeze: Calculates the amount of refreezing undergone in the subsurface
%column after the meltwater has percolated.
%   syntax:
%       [psnowc, psnic, pslwc, ptsoil, zrfrz, pts ]...
%           = refreeze (psnowc, psnic, pslwc, ptsoil, c)
%   input:
%         psnowc, psnic, pslwc - vectors of respectively snow, ice and
%         liquid water part for each subsurface layer (mm weq). Their sum
%         for each layer should always be equal to the layer fixed water
%         eq. thickness.
%
%         ptsoil - vector of subsurface temperature (K)
%
%         c - Structure containing all the physical, site-dependant or user
%         defined variables.
%
%   output:
%          updated prhofirn, ptsoil, psnowc, psnic, pslwc
%          zrfrz - amount of water (mm eq) refrozen during the time step
%          pts - surface temperature
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to IDL by
%   Robert S. Fausto and later to Matlab by Baptiste Vandecrux
%   (bava@byg.dtu.dk).
% ========================================================================

%======================================================
%  Here we do the refreezing based on the cold content
% of each layer converting mass from liquid to ice.
%======================================================
zrfrz=zeros(c.jpgrnd,1);
zpotref=zeros(c.jpgrnd,1);
coldcontent=zeros(c.jpgrnd,1);
cfrozen=zeros(c.jpgrnd,1);

% Determine layer cold content and convert mass from liq to ice
cfrozen = psnowc + psnic;
coldcontent = max(0, (c.T_0-ptsoil));

% update summer 2013 by BV
% only a fraction C_k * dt of the potential refrozen water is allowed to
% refreeze during the time step dt. See Illangaskare et al. 1990, eq. 13
zpotref = c.Ck * ...
    coldcontent .* cpiceF(ptsoil) .* cfrozen / c.L_fus;
zrfrz= min(zpotref , pslwc);  
    
pslwc =  pslwc - zrfrz;
psnic = psnic + zrfrz;

% Update temperature with c.latent heat of freezing
ptsoil = ptsoil + ...
    zrfrz.* c.L_fus./( cpiceF(273.15).*(psnic-zrfrz +...
    psnowc) );
% pts = ptsoil(1);
end
