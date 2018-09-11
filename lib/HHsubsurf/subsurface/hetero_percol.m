function [pslwc] = hetero_percol (prhofirn, psnowc , psnic, pslwc, pdgrain, c)
% hetero_percol
%   syntax:
% [prhofirn, psnowc , psnic, pslwc, pdgrain] =...
%         hetero_percol (prhofirn, psnowc , psnic, pslwc, ...
%         pdgrain, c)
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
%   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk)
%=========================================================================
avail_water = zeros(size(prhofirn));

for jk = 1:c.jpgrnd-1;
    if ( c.calc_CLliq )
         liqmaxloc = CLliqF(prhofirn(jk),c);
      else
         liqmaxloc = c.liqmax;
    end

    liqmaxM = liqmaxloc*c.rho_water/c.rho_ice*(c.rho_ice/prhofirn(jk) - 1);
    potret    = max( liqmaxM* psnowc(jk) , 0 );
    liqexcess = pslwc(jk) - potret;
    avail_water(jk) = max ( liqexcess , 0 );
    
    %if there is available water
    if avail_water(jk)>c.smallno
        % we draw from a binomial distribution of mode hetero_percol_p to
        % know whether this layer leads to heterogeneous percolation
        if binornd(1,c.hetero_percol_p)
            fprintf('Piping occuring\n')
            percol_water = c.hetero_percol_frac * avail_water(jk);
            %note: dflux is positive when it leaves the layer
            
            %determine random depth between current layer and maximum depth
            %range
%             thickness_weq = psnowc + snic +slwc;
            thickness_act = psnowc.*(c.rho_water./prhofirn)+ psnic.*(c.rho_water/c.rho_ice);
            depth_act=zeros(size(thickness_act));
            depth_act(1)= thickness_act(1)/2 ;
            for i = 2:size(thickness_act,1)
                depth_act(i)=sum(thickness_act(1:i-1)) + thickness_act(i)/2 ;
            end

            depth_current = depth_act(jk);
            depth_dest = c.hetero_percol_dist * rand(1,1) + depth_current;
            % find index of layer at closest depth to the destination of
            % percolation
            [~, ind_dest] = min(abs(depth_act-depth_dest));
            
            fprintf('from %0.2f m deep to %0.2f m\n',depth_current,depth_dest)

            for ii = jk:ind_dest-1
                %here we test all the layer through which the pipe travels
                if ThetaF(pslwc(ii+1), psnowc(ii+1), prhofirn(ii+1), c) >= 1
                    %if there is a saturated layer piping can't go through
                    %it
                    fprintf('sat - ')
                    break
                elseif psnic(ii+1)> c.ice_threshold
                    % plus relative to frozen mass 
                    %if there is too much ie in the next layer
                    fprintf('ice - ')
                    break
% PERMEABILITY? h?
%                 elseif pdgrain(ii)-pdgrain(ii+1) < c.crit_diff_grain
%                     fprintf('grain - ')
%                     break
                end
            end
            if ind_dest ~= ii
                fprintf('Piping stopped at layer %i instead of \n',ii, ind_dest)
            end
            ind_dest = ii;
            
            jj = ind_dest;
            while and(jj>jk , percol_water > c.smallno)
                % Calculate water in destination layer, when this is at saturation (Theta = 1):
                plsat = psnowc(jj)*c.rho_water/c.rho_ice *...
                    (c.rho_ice/prhofirn(jj) - 1);
                % Do not allow flux to be greater than plsat-pl in next layer.
                dflux = min (max(0, plsat-pslwc(jj)), percol_water);
                                
                pslwc(jk) = pslwc(jk) - dflux;
                pslwc(jj) = pslwc(jj) + dflux;
                percol_water = percol_water - dflux;
                jj = jj-1;
            end
        end
    end
end
end
