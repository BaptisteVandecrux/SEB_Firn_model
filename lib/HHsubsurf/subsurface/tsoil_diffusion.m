function [ptsoil] = tsoil_diffusion (pts, pgrndc, pgrndd, ptsoil, c)
%   tsoil_diffusion: Update subsurface temperatures ptsoil based on previous 
%   time step coefficients for heat conduction
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================

%         ! Upper layer
ptsoil(1)= pts;

%         ! Lower layers
% % update BV2017
% for jk = 2:c.jpgrnd
%     ptsoil(jk) = pgrndc(jk) + pgrndd(jk) * ptsoil(jk);
% end
% WHY??
% original 
for jk = 1:c.jpgrnd-1
    ptsoil(jk+1) = pgrndc(jk) + pgrndd(jk) * ptsoil(jk);
end
end
