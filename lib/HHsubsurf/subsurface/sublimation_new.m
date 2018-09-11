
function [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]...
    = sublimation_new (zsn, psnowc, psnic, pslwc, pdgrain ...
    , prhofirn, ptsoil, snowbkt, c)
% sublimation_new: Routine that removes sublimation from the first layer of
% the column. Update BV2017: The proportion of ice and snow that is 
% sublimated is given by the proportion of ice
% and snow present in the layer. If the first layer
% becomes smaller than lim_old_lay, it is merged with the second
% layer,another layer is splitted elsewhere and sublimation continues on 
% this new merged layer. Material from the snowbkt is taken before anything
% else.
%	
%   input:
%         zsn - mass (m weq) of material that needs to be sublimated
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
%   This script was developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
%   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
%=========================================================================

if (zsn < 0)  % zsn negative means net upward snowfall+sublimation
    zsnout = abs(zsn) ; % zsnout is the mass the leaves column upwards
    
    % first  we sublimate the content of the snow bucket
    if snowbkt > c.smallno
        pot_subl = min(snowbkt, zsnout);
        snowbkt = snowbkt - pot_subl;
        zsnout = zsnout - pot_subl;
    end
    
    %then comes the upper  layers of the column
    while (zsnout > c.smallno )

        % Update BV2017: Snow and ice are now melted simultaneously
        zdel = min(zsnout, psnowc(1) + psnic(1));
        snow_melt = psnowc(1)/(psnowc(1)+psnic(1)) * zdel;
        ice_melt = psnic(1)/(psnowc(1)+psnic(1)) * zdel;

        psnowc(1) = psnowc(1) - snow_melt;
        psnic(1) = psnic(1) - ice_melt;
        zsnout =zsnout - zdel;

        if psnowc(1) + psnic(1) + pslwc(1) < c.lim_new_lay
            if (psnowc(2) + psnowc(1))> c.smallno
                %then merging layer 1 and 2 and placing it in layer 2
                pdgrain(2) = (psnowc(1) * pdgrain(1) + psnowc(2) * pdgrain(2)) ...
                   /(psnowc(2) + psnowc(1));
                prhofirn(2) = (psnowc(2) + psnowc(1)) ...
                   / (psnowc(1)/prhofirn(1) +  psnowc(2)/prhofirn(2));
%                Update PLA
                ptsoil(2) = ((psnic(1) + psnowc(1) + pslwc(1)) * ptsoil(1) ...
                    + (psnic(2) + psnowc(2) + pslwc(2)) * ptsoil(2)) ...
                    /(psnic(1) + psnowc(1) + pslwc(1) + psnic(2) + psnowc(2) + pslwc(2));


                psnowc(2) = psnowc(2) + psnowc(1);
                psnic(2) = psnic(2) + psnic(1);
                pslwc(2) = pslwc(2) + pslwc(1);
            else
                pdgrain(2) = -999;
                prhofirn(2) = -999;
            end
            
            % now top layer is free and another layer can be splitted
            % elsewhere
            [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = ...
                split_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil, c);
        end
    end
end
end