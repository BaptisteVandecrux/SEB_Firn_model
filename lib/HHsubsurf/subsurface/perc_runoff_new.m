function [prhofirn, psnowc , psnic, pslwc, pdgrain, zrogl] =...
    perc_runoff_new (prhofirn, psnowc , psnic, pslwc, ...
        pdgrain, c)

% perc_runoff_new: Calculates meltwater percolation and runoff in the column
% either according to a standard bucket scheme or to a Darcy-like bucket
% scheme. Update BV2017: no mass shift anymore, the water is
% moving from one layer to the next and layer's total mass is allowed to change.
% At the end of the routine, it is checked that no layer has a mass lower
% than lim_old_lay. If one is found too small, it is merged with its
% underlying layer and another layer is split elsewhere.
%
%   syntax:
%   [prhofirn, psnowc , psnic, pslwc, ptsoil , pdgrain, pTdeep, zrogl] =...
%     perc_runoff_new (prhofirn, psnowc , psnic, pslwc, ptsoil , ...
%           pdgrain, pTdeep, c.ElevGrad, zrogl, c)
%	
%   input:
%         prhofirn - vector of firn (snow) density (kg/m^3) 
%
%         psnowc, psnic, pslwc - vectors of respectively snow, ice and
%         liquid water part for each subsurface layer (m weq).
%
%         ptsoil - vector of subsurface temperature (K)
%
%         pdgrain - Vector of layer grain size. see graingrowth function.
%
%         pTdeep - lower boundary temperature (K). Defined by the constant
%         T_ice_AWS so far.
%
%         zrogl - total amount of run off (mm weq)
%
%         c - Structure containing all the physical, site-dependant or user
%         defined variables.
%
%   output:
%          updated [prhofirn, psnowc , psnic, pslwc, ptsoil , pdgrain, zrogl]
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================

%======================================================
%Here we do liquid water percolation and runoff
%======================================================
% *** Calculate potential Darcy fluxes for all columns ***
%  (if do-no-darcy, then fluxes will be calcuc.lated on-the-fly in loops below)
% This calculation makes sure potential darcy flux across interfaces does
% not exceed the water amount in the upper of the two layers

if (~c.do_no_darcy)
    [darcy_fluxes] = calc_darcy_fluxes (pslwc, psnowc , psnic, pdgrain ...
        , prhofirn, c);
end

%Update BV2017: t_runoff is now calculated outside of the subsurface scheme

% Bottom layer: Receive from above. Give excess liquid to runoff.
jk = c.jpgrnd;
zrogl = 0;
    % Remove runoff from layer (and add to runoff-box zrogl)
    % Adjust bottom layer interface due to removal of mass
    % PLA Darcy 2016
     if ( c.calc_CLliq )
             liqmaxloc = CLliqF(prhofirn(jk),c);
          else
             liqmaxloc = c.liqmax;
     end

    liqmaxM = liqmaxloc*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk) - 1);
    potret    = max( liqmaxM* psnowc(jk) , 0 );
    % Calculate liqexcess from layer water. Make sure it is not negative:
    liqexcess = max ( pslwc(jk) - potret , 0 );
    liqro = liqexcess/c.t_runoff * c.zdtime;
    % Take runoff from water content and give to runoff box (Update PLA)
    zrogl    = zrogl    + liqro  ;
    pslwc(jk) = pslwc(jk) - liqro;

for jk = c.jpgrnd-1:-1:1
        % BV2017 removing percolation blocking criteria
        if ThetaF(pslwc(jk+1), psnowc(jk+1), prhofirn(jk+1), c) >= 1
            % if next layer already saturated, then we do not allow flux 
            % into that layer, but instead activate runoff
            dflux = 0;
            do_runoff = 1 ;
        else
            if ( c.do_no_darcy )
                % Potential Darcy-like flux in case of do-no-Darcy,
                % basically just all liqexcess:
                % PLA Darcy 2016
                 if ( c.calc_CLliq )
                         liqmaxloc = CLliqF(prhofirn(jk),c);
                      else
                         liqmaxloc = c.liqmax;
                 end

                liqmaxM = liqmaxloc*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk) - 1);
                potret    = max( liqmaxM* psnowc(jk) , 0 );
                liqexcess = pslwc(jk) - potret;
                darcy_fluxes(jk) = max(liqexcess , 0 );
            end

            % Calculate water in next layer, when this is at saturation (Theta = 1):
            plsat = psnowc(jk+1)*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk+1) - 1);

            % Do not allow flux to be greater than plsat-pl in next layer.
            % Also, do not allow it to be negative
            dflux = max ( min( darcy_fluxes(jk) , plsat-pslwc(jk+1) ) , 0  );

                % Update BV2017: if no snow to hold liquid water then runoff is
                % turned on
            if or( darcy_fluxes(jk) >= plsat-pslwc(jk+1),...
                    psnowc(jk)<c.smallno)
                % There is enough darcy flow to fill up layer below (and perhaps more)
                % so runoff should be activated (on what's left here after
                % we have moved down enough to fill layer below):
                do_runoff = 1;
            else
                % No runoff from this layer, since darcy flow is not large enough
                % to fill up layer below
                do_runoff = 0;
            end
        end

        if ( dflux > 0 )
            % Yes: Darcy flow to next layer
               % Update BV2017: Now temperature is not updated. water
               % percolating is at 0degC cannot change temperature until phase
               % change occurs in refreezing scheme
            pslwc(jk+1) = pslwc(jk+1) + dflux;
            pslwc(jk) = pslwc(jk) - dflux ;
        end
        
        if and(do_runoff, ~c.avoid_runoff)
            % Yes: Remove runoff from layer (and add to runoff-box zrogl)
            % PLA Darcy 2016
             if ( c.calc_CLliq )
                     liqmaxloc = CLliqF(prhofirn(jk),c);
                  else
                     liqmaxloc = c.liqmax;
             end

            liqmaxM = liqmaxloc*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk) - 1);
            potret    = max( liqmaxM* psnowc(jk) , 0);
            
            % Calculate liqexcess from layer water. Make sure it is not negative:
                % Update BV2017: since dflux has already left slwc(jk) so no
                % need to remove it anymore (compared to old version)
            liqexcess = max ( pslwc(jk) - potret , 0 );
                       
                %Update BV2017: for all layer, if there is no snow to hold
                %liquid water then runoff is instantaneous ( not allowing
                %subglacial lakes to form
            if (psnowc(jk) < c.smallno)
                % Here in layer: If there is no snow, run off immediately
                % in other word surface runoff is instantaneous
                liqro = liqexcess;
            else
                Theta = ThetaF( pslwc(jk),psnowc(jk),prhofirn(jk), c);
%                 liqro_darcy = kF( Theta,pdgrain(jk),prhofirn(jk),psnic(jk) ,psnowc(jk), c) * c.ElevGrad; 
                % old version based on Zuo and Oerlemans (1996)
                liqro_darcy  = liqexcess/c.t_runoff * c.zdtime;
                
                liqro = min(liqro_darcy,liqexcess);
                pore_space = psnowc(jk)*c.rho_water*(1/prhofirn(jk)-1/c.rho_ice);

                %Update PLA
                if pslwc(jk) - liqro > pore_space
                    liqro = pslwc(jk) - pore_space;
                end
                    
%                 fprintf('Runoff %f %f %f\n',liqro_darcy, liqro_zo, liqro);
            end
            
            % Take runoff from water content and give to runoff box
            zrogl    = zrogl    + liqro  ;
            pslwc(jk) = pslwc(jk) - liqro;
        end
end
end
