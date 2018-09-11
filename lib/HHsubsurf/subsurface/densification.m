function [prhofirn, dH_comp, dV] = densification (pslwc, psnowc , prhofirn, ptsoil, c)
%   densification: updates the subsurface firn/snow density with
%   gravity-driven densification.
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================

%   ! PLA densification (Feb 2015)
%   ! Densification of layers (RSF+PLA Jan-Feb 2015):
%   ! Vionnet et al. (GMD, 2012) densification eqs 5-9

    %            ! Pressure in Pa (evaluated at layer mid-point).
    %            ! We have neglected surface slope in eqn 6.

    
    sigma_pres(1:c.jpgrnd,1) = c.cmid * c.rho_water * c.g ;
    
    f1f2(1:c.jpgrnd,1) = 4./( 1 + 60.*(pslwc./ (psnowc + c.smallno)) ...
        .* (prhofirn ./ c.rho_water) ); % ! f2 = 4
    
    eta_firn(1:c.jpgrnd,1)   = f1f2 .* 7.62237e+6 .* (c.a_dens * prhofirn / 250 ) ...
        .* exp( 0.1.*(c.T_0 - ptsoil) + c.b_dens * 0.023.*prhofirn );
    prhofirn_old = prhofirn;
    prhofirn = prhofirn + c.zdtime*prhofirn .* sigma_pres ./eta_firn;

    % Update Baptiste 2018
%     MO = NaN(size(prhofirn));
%     MO(prhofirn <= 550)  = max(0.25, 1.042 - 0.0916 * log(c.accum_AWS));
%     MO(prhofirn > 550)  = max(0.25, 1.734 - 0.2039 * log(c.accum_AWS));
% 
%     C = NaN(size(prhofirn));
%     C(prhofirn <= 550)  = 0.07;
%     C(prhofirn > 550)  = 0.03;
%     
%     Ec = 60*1000; % Jmol?1
%     Eg = 42.4*1000; % J mol?1
%     prhofirn = prhofirn + c.zdtime.*...
%         MO .* C .* c.accum_AWS * c.rho_water/365/24/3600 .* c.g .* (c.rho_ice - prhofirn) ...
%         .* exp(-Ec./c.R./ptsoil + Eg./c.R./ptsoil);

% Arthern
%     c.R = 8.314;
%     Ec = 60*1000;
%     Eg = 42.4*1000;
%     m0 = 1; %1.042 - 0.0916*log(acc);
%     m1 = 1; %1.734 - 0.2039*log(acc);
%     % since the standard time unit in the code is second
%     % we need to convert accumulation rates to kg/s
%     c_0 = 0.07 * c.mean_accum*c.rho_water/365/24/3600*c.g .* exp( - Ec / c.R ./ptsoil + Eg / c.R / c.Tdeep);
%     c_1 = 0.03 * c.mean_accum*c.rho_water/365/24/3600*c.g .* exp( - Ec / c.R ./ptsoil + Eg / c.R / c.Tdeep);
%     
%     ind1 = ( prhofirn <= 550);
%     prhofirn(ind1) = prhofirn(ind1) + c.zdtime * m0.* c_0(ind1).*(c.rho_ice - prhofirn(ind1));
%     ind2 = prhofirn > 550;
%     prhofirn(ind2) = prhofirn(ind2)  + c.zdtime * c_1(ind2) .* m1 .*(c.rho_ice - prhofirn(ind2));
    
    prhofirn(prhofirn > c.rho_ice)= c.rho_ice;

    dV = psnowc * c.rho_water ./prhofirn_old - psnowc * c.rho_water ./prhofirn;
    dH_comp = sum(dV); %in real m, here positive when compaction o
end