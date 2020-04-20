function [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = ...
            split_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil, c)
% split_layer: This function finds the layer that is the best to be
% splitted. In most cases, since gradients in temperature, water content
% and density are greater close to the surface, it will select the
% uppermost layer that has enough mass to be splitted (mass >
% lim_new_lay*2). If all layers in the column are too small to be split,
% then the model will create a new layer at the bottom of the column with
% 10 m of ice and a user defined temperature of Tdeep.
% WARNING: the top layer needs to be free: its content will be discarded
% during the renumbering of the column.
%	
%   input:
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
%         c - Structure containing all the physical, site-dependant or user
%         defined variables.
%
%   output:
%          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil]
%
%   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
%   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
%=========================================================================   

    thickness_weq = psnic + psnowc + pslwc;
    
    if and(sum(thickness_weq > c.lim_new_lay*2) > 0,...
            sum(thickness_weq) > c.min_tot_thick)
%         if there are some layers that have are thick enough to be
%         splitted into two layers of which thickness would be higher than
%         the limit thickness for new layer c.lim_new_lay

        depth_weq = cumsum(thickness_weq);
        mid_point_weq = depth_weq - thickness_weq./2;

        %first criterion: Layer that are thicker are more likely to be splitted
        crit_1 = interp1([0 c.lim_new_lay  c.thick_crit], [0 0 1], thickness_weq,'linear', 1);

        %second criterion: Shallow layers are more likely to be splitted    
        crit_2 = interp1([0 c.deep_firn_lim], [1 0], mid_point_weq,'linear', 0);

        crit = (0.6*crit_1 + 0.4*crit_2) / 2;
        
        %all layer that are too thin to be splitted are taken out
        too_thin = thickness_weq./2 < c.lim_new_lay;
        crit(too_thin) = -1;
        
        crit(1) = 0; %We don't want to split the first layer anyway
        [~, i_split] = max(crit);
%         fprintf('Splitting layer %i\n',i_split);

        psnowc(1:i_split-1) = psnowc(2:i_split);
        psnic(1:i_split-1) = psnic(2:i_split);
        pslwc(1:i_split-1) = pslwc(2:i_split);
        pdgrain(1:i_split-1) = pdgrain(2:i_split);
        prhofirn(1:i_split-1) = prhofirn(2:i_split);
        ptsoil(1:i_split-1) = ptsoil(2:i_split);
        % now layer i_split and i_split-1 are the same
        
        thickness_new = min(c.max_lay_thick, thickness_weq(i_split)/2);
        thickness_left = thickness_weq(i_split) - thickness_new;
        r_new = thickness_new / thickness_weq(i_split);
        r_left = thickness_left / thickness_weq(i_split);

        psnowc(i_split-1) = r_new*psnowc(i_split);
        psnic(i_split-1) = r_new*psnic(i_split);
        pslwc(i_split-1) = r_new*pslwc(i_split);
        % prhofirn, ptsoil and pdgrain were already identical

        psnowc(i_split) = r_left*psnowc(i_split);
        psnic(i_split) = r_left*psnic(i_split);
        pslwc(i_split) = r_left*pslwc(i_split);
        % now layer is split
    else
%       if all layers are too small, then ice is taken from below.
        psnowc(1:end-1) = psnowc(2:end);
        psnic(1:end-1) = psnic(2:end);
        pslwc(1:end-1) = pslwc(2:end);
        pdgrain(1:end-1) = pdgrain(2:end);
        prhofirn(1:end-1) = prhofirn(2:end);
        ptsoil(1:end-1) = ptsoil(2:end);
        
        psnowc(end) = 0;
        psnic(end) = c.new_bottom_lay;
        pslwc(end) = 0;
        pdgrain(end) = 1;
        prhofirn(end) = 1;
        ptsoil(end) = c.Tdeep;
%         disp('ADDING NEW LAYER')
    end
end