function [pdgrain] =  graingrowth ( pslwc, psnowc , pdgrain, c)
%   graingrowth: Update the subsurface grain size vector after grain growth
%   due to metamorphism.

%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================

%   ! PLA Darcy (may 2016)
pdgrain( 1:c.jpgrnd) = pdgrain( 1:c.jpgrnd) + ...
    c.zdtime*dgraindtF(pdgrain( 1:c.jpgrnd),...
    pslwc( 1:c.jpgrnd),psnowc(1:c.jpgrnd), c);
end