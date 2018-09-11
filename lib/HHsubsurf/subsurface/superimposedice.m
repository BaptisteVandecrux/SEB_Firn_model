function [ptsoil , psnic, pslwc, zsupimp]...
    =  superimposedice (prhofirn, ptsoil            ...
    , psnowc  , psnic, pslwc, zso_cond, c )
%superimposedice: Calculates the amount (zsupimp in mm eq) of liquid water 
% that is refrozen on top of perched ice layers. New version from 2016. 
% Does not use slush bucket.
%   syntax:
% [prhofirn, ptsoil , psnowc  , psnic, pslwc, zsupimp]...
%     =  superimposedice (prhofirn, ptsoil            ...
%     , psnowc  , psnic, pslwc, zso_cond, c )
%	
%   input:
%         prhofirn - vector of firn (snow) density (kg/m^3) 
%
%         psnowc, psnic, pslwc - vectors of respectively snow, ice and
%         liquid water part for each subsurface layer (mm weq). Their sum
%         for each layer should always be equal to the layer fixed water
%         eq. thickness.
%
%         ptsoil - vector of subsurface temperature (K)
%
%         zso_cond - Subsurface thermal conductivity. See ice_heats
%         function.
%
%         c - Structure containing all the physical, site-dependant or user
%         defined variables.
%
%   output:
%          updated [psnowc , psnic, pslwc, ptsoil]
%          zsupimp - total amount of formed superimposed ice (mm weq)
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================

% --------------------------------------
% ----- SUPERIMPOSED ICE FORMATION -----
% --------------------------------------

% In the "Darcy-version", superimposed ice formation is changed in these respects:
%   1 Water used is just slwc (and not a slush bucket)
%   2 This simplifies things because no mass is added, subtracted or shifted.
%     Just converted from slwc to snic.
%   3 The new ice used to be added to the cold ice layer below. Now it is simply
%     converted to ice in the layer, where the water resides (giving energy to
%     the cold layer below).
%   4 SIF can oc.ccur several places in the column, as there is no c.longer a slush
%     bucket. slwc is now allowed to be larger than potret (the excess just runs off
%     with the Zuo&Oerlemans timescale as the sluch water used to) and it is this
%     stock of slwc that is used for SIF.

% No SIF allowed onto inifinte sublayer, so only
% search for suitable water layer in 1,c.jpgrnd-1 (ice layer below is in 2,c.jpgrnd):
% jk is the layer our water is in. We then look at jk+1 to see if
% we can superimpose:
% PLA-BV 2016: volume (instead of mass) weighted average of density
rho_next(1:c.jpgrnd-1) = (psnic(2:c.jpgrnd)+psnowc(2:c.jpgrnd))./...
    ( psnic(2:c.jpgrnd)/c.rho_ice + psnowc(2:c.jpgrnd) ./ prhofirn(2:c.jpgrnd) );

nextisfrozen(1:c.jpgrnd-1) = (rho_next(1:c.jpgrnd-1) >= c.rho_pco );
isfrozen(1)=0;
isfrozen(2:c.jpgrnd) = nextisfrozen(1:c.jpgrnd-1);
% The next layer has bulk density high enough that we will consider it as
% an ice layer. Therefore we calcuc.late SI-formation:

snowV = 0.5*psnowc(isfrozen==1)  .*c.rho_water./prhofirn(isfrozen==1);
iceV  = 0.5*psnic(isfrozen==1)  *c.rh2oice;
totalV = snowV + iceV;
zx1 = snowV ./ totalV;
zx2  = iceV   ./ totalV;
% update BV2017 zsn_cond is now vector
ki= 1./( zx1./zsn_condF(prhofirn(isfrozen==1), c) + zx2./zso_cond(isfrozen==1) );
dTdz =  (c.T_0-ptsoil(isfrozen==1))  ./ totalV ;

% The potential SIformation (Only interested in positive SIF):
potSIform = max(0, ki.*dTdz / (c.rho_water*c.L_fus) * c.zdtime) ;

% Make as much SI as there is water available in the layer above the ice (jk)
SIform = min(pslwc(nextisfrozen==1),potSIform);
zsupimp = zeros(size(pslwc));
zsupimp(nextisfrozen==1) = SIform;

% Convert water to ice
pslwc(nextisfrozen==1) = pslwc(nextisfrozen==1) - SIform;
%BV 2016 corrected to avoid underflow
psnic(nextisfrozen==1) = psnic(nextisfrozen==1) + SIform;
 
% Update temperature of ice layer
ptsoil(isfrozen==1) = ptsoil(isfrozen==1) + ...
    SIform * c.L_fus ./ (c.cdel(isfrozen==1).*cpiceF(ptsoil(isfrozen==1)));
end

