function [dgraindtF] = dgraindtF(dg,pl,ps, c)
%     ! The F in the name stands for "function"
%     ! Calculates d(grainsize-diameter) /dt (in mm/s) as in Hirashima (2010), using saturated
%     ! and unsaturated rates from Brun (1989), Tusima (1978) and Katsuhima (2009).

%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================
%     ! Mass liquid water content in %
L = pl./(ps+c.smallno).*100;

%     ! Brun derives a function for grain growth that increases as L^3 (and divides by d^2).
%     ! Katushima uses this for L<=10%. Beyond 10 % Katushima uses the minimum value of the
%     ! Brun formuc.lation and the constant-divided-by-d^2 found by Tusima.
%     ! However, the Tusima-constant is lower than the L^3 function at 10% so that would lead
%     ! to a discontinuous drop when going beyond 10%.
%     ! Instead, we do as we also think Hirashima does "Therefore, to consider grain growth
%     ! in the water-saturated layer, the grain coarsening model of Tusima (1978) was used
%     ! as an upper boundary to saturated grain growth following Katsushima et al. (2009a)":

%     ! Take the L^3 function by Brun up until this becomes larger than the Tusima-constant:
Brun   = 2./pi.* ( 1.28e-8 + 4.22e-10 .* L.^3);
Tusima = 6.94e-8;
%     ! d (grainsize-diamter) /dt in mm/s
dgraindtF = 1./(dg.^2).*min(Brun,Tusima);
end
