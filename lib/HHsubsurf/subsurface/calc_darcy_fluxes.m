function [darcy_fluxes] = calc_darcy_fluxes (pslwc, psnowc , psnic,...
    pdgrain , prhofirn, c)

%calc_darcy_fluxes: Calculates the amount of water (mm weq) that each layer
%can transmit to the next one according to a Darcy flow.
%
%   syntax:
%   [darcy_fluxes] = calc_darcy_fluxes (pslwc, psnowc , psnic,...
%            pdgrain , prhofirn, c)
%	
%   input:
%         psnowc, psnic, pslwc - vectors of respectively snow, ice and
%         liquid water part for each subsurface layer (mm weq). Their sum
%         for each layer should always be equal to the layer fixed water
%         eq. thickness.
%
%         pdgrain - Vector of layer grain size. see graingrowth function.
%
%         prhofirn - vector of firn (snow) density (kg/m^3) 
%
%         c - Structure containing all the physical, site-dependant or user
%         defined variables.
%
%   output:
%         darcy_fluxes - vector containing the amount of water (mm weq) that each layer
%         can transmit to the next one according to a Darcy flow.
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%
%   Aug. 2016
%=========================================================================

darcy_fluxes = zeros(c.jpgrnd,1);
for jk = 1:c.jpgrnd-1
    % Special treatment of top layer in case of no snow
    if (jk == 1 && (psnowc(jk) < c.smallno))
        % In case of no snow, the calculation below of K crashes. Instead we say
        % * if layer below may receive water, it gets all there is in layer 1
        % * if layer may not receive water, it all runs off.
        %
        % As long as we set darcy_fluxes(1) = pslwc(1) then this will
        % be taken care of in the perc_runoff code.
        darcy_fluxes(jk) = pslwc(jk);
    else
        % Other layers do (should!) not come into a situation where there is only ice and water
        qlim = qlimF(                      ...
            pslwc(jk)   , pslwc(jk+1)    , ...
            psnowc(jk)  , psnowc(jk+1)   , ...
            psnic(jk)   , psnic(jk+1)    , ...
            prhofirn(jk), prhofirn(jk+1) , ...
            pdgrain(jk) , pdgrain(jk+1), c );
        Theta1 = ThetaF( pslwc(jk)  ,psnowc(jk)   ,prhofirn(jk), c );
        Theta2 = ThetaF( pslwc(jk+1),psnowc(jk+1) ,prhofirn(jk+1), c);
        h1     = hHirF( Theta1, pdgrain(jk), c);
        h2     = hHirF( Theta2, pdgrain(jk+1), c);

        delz1 = c.rho_water*(psnowc(jk)  /prhofirn(jk)   + psnic(jk)  /c.rho_ice); % layer 1
        delz2 = c.rho_water*(psnowc(jk+1)/prhofirn(jk+1) + psnic(jk+1)/c.rho_ice); % layer 2
        delz = (delz1+delz2)/2;             % Midpoint-to-midpoint distance
        
        dhdz = (h2-h1)/delz;
        
        k1 = kF( Theta1,pdgrain(jk)  ,prhofirn(jk)  ,psnic(jk)  ,psnowc(jk), c);
        k2 = kF( Theta2,pdgrain(jk+1),prhofirn(jk+1),psnic(jk+1),psnowc(jk+1), c);
        k = (k1+k2)/2;                      % Average k between the two layers

        % Initial flux according to Hirashima eqn (1)
        q0 = max(k*(dhdz + 1),0);          %
        
        % Total time-step flux according to Hirashima eqn (23)
        % and make sure it doesn't exceed available water in
        % the upper of the two layers:
        if (qlim > 0)        
            darcy_fluxes(jk) = min( pslwc(jk) , qlim*(1-exp(-q0/qlim*c.zdtime)) );
        else
            darcy_fluxes(jk) = 0;
        end
        
    end
end

end
