function [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]...
    = snowfall_new (zsn, psnowc, psnic, pslwc, pdgrain ...
    , prhofirn, ptsoil, pts, snowbkt, zraind, zsnmel, c)
% snowfall_new: Routine that adds new fallen snow (net from 
% sublimation)on top of the subsurface column. It is first accumulated in 
% snowbkt and only when snowbkt exceed a threshold then a new layer is
% created. Update BV2017: No mass shift anymore. When a new layer is created, two 
% layers at depth are merged.
%	
%   input:
%         zsn - amount (m weq) of new fallen snow net of sublimationaning
%         that it can be negative if sublimation is greater than
%         snow accumulation).
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
%         snowbkt - Fresh snow bucket where new snow is added. When its
%         content is greater than lim_new_lay, then a new layer is created.
%
%         c - Structure containing all the physical, site-dependant or user
%         defined variables.
%
%   output:
%          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]
%
%   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
%   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
%=========================================================================

if (zsn > 0) % ! zsn means snowfall in timestep

    %snowfall added to the fresh snow bucket
    snowbkt = snowbkt + zsn;
    
    if snowbkt > c.lim_new_lay
        %enough material for a new layer           
%         if sum(psnic+psnowc+pslwc)> c.lim_max_column
%         % Update BV2017
%         % if the column (excluding last layer) is larger than a threshold, 
%         % then last layer is discarded, all the layers are shifted downward
%         % and the new snow layer is put in first layer.
%             psnowc(2:end) = psnowc(1:end-1);
%             psnic(2:end) = psnic(1:end-1);
%             pslwc(2:end) = pslwc(1:end-1);
%             pdgrain(2:end) = pdgrain(1:end-1);
%             prhofirn(2:end) = prhofirn(1:end-1);
%             ptsoil(2:end) = ptsoil(1:end-1);
%         else
            [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = ...
                merge_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil,c);
%         end
%         disp('SNOWFALL')
%         disp(snowbkt)
        psnowc(1) = snowbkt;
        psnic(1) = 0;
        pslwc(1) = 0;

        ptsoil(1) = min(pts,c.T_0);
        prhofirn(1) = c.rho_fresh_snow;
        pdgrain(1) = c.dgrainNew;
        snowbkt = 0;
    end
end

% Update BV 2017
% if there is rain or melt, the fresh snow bucket is added to the snow
% compartment of the first layer
if and(or(zraind>c.smallno,zsnmel>c.smallno),snowbkt>c.smallno)
    % mass-weighted average for temperature
    ptsoil(1) = (ptsoil(1)*(psnowc(1)+pslwc(1)+psnic(1)) + min(pts,c.T_0)*snowbkt)...
        /(psnowc(1)+pslwc(1)+psnic(1)+ snowbkt);

    %snow-mass-weighted average for grain size
    pdgrain(1) = (psnowc(1) * pdgrain(1) + snowbkt * c.dgrainNew) ...
           /(psnowc(1) + snowbkt);
    
    % volume-weighted average for density
    prhofirn(1) = (psnowc(1) + snowbkt) ...
       / (snowbkt/c.rho_fresh_snow +  psnowc(1)/prhofirn(1));
 
   psnowc(1) = psnowc(1) + snowbkt;
   snowbkt = 0;       
end
    
end