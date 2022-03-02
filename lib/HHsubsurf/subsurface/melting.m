function [psnowc, psnic, pslwc, ptsoil, snowbkt] = ...
    melting(psnowc, psnic, pslwc, zsnmel, snowbkt, ptsoil, prhofirn, c)
% melting: Transform an amount zsnmel (m eq) of frozen material into liquid
% water. Starts at the surface and continues downward as long as more
% material needs to be melted. Update BV2017: The proportion of ice and snow  
% that is melted is given by the proportion of ice and snow present in the 
% layer. Material from the snowbkt is taken before anything else. The
% volume of what is being melted is calculated.
%	
%   input:
%         zsnmel - mass (m weq) of material that needs to be melted
%         prhofirn - vector of firn (snow) density (kg/m^3) 
%         psnowc, psnic, pslwc - vectors of respectively snow, ice and
%         liquid water part for each subsurface layer (m weq).
%         ptsoil - vector of subsurface temperature (K)
%         snowbkt - Fresh snow bucket where new snow is added. When its
%         content is greater than lim_new_lay, then a new layer is created.
%         c - Structure containing all the physical, site-dependant or user
%         defined variables.
%
%   output:
%          updated [psnowc, psnic, pslwc, pdgrain, prhofirn, ptsoil, snowbkt]
%           vol_melted - volume of the material that was melted
%
%   This script was originally developped by Baptiste Vandecrux (bava@byg.dtu.dk), 
%   Peter L. Langen (pla@dmi.dk) and Robert S. Fausto (rsf@geus.dk).
%=========================================================================

% First we handle melting. Make sure all melt energy is used. Start by melting in
% top layer. If necessary go to deeper layers.
zsnout = zsnmel;

jk = 0;

% first  we melt the content of the snow bucket
if (snowbkt > c.smallno && zsnout > c.smallno)
    zdel = min(snowbkt, zsnout);
    snowbkt = snowbkt - zdel;
    zsnout = zsnout - zdel;
    pslwc(1) = pslwc(1) + zdel;
end

%         ! Enter if there is any melting
while (zsnout > c.smallno)
    %  Start at top and work downward if necessary
    jk = jk + 1 ;
    
    % Update BV2017
    % Normally, if melt is prescribed at the surface, it implies that the
    % top layer is at melting point. Since we allow the melting of multiple
    % layer when needed, it can occur that we need to melt a subfreezing layer.
    % In that situation, the melting energy is used to warm that layer.
    % This warming occurs without phase change.
    
    % How much energy is needed to bring the layer to melting point
    deltaT = c.T_0 - ptsoil(jk);
    % volumetric heat capacity of the layer in J/m3/K
    heat_capa_vol = zsn_capaF(prhofirn(jk),ptsoil(jk));
    volume = (psnic(jk) + psnowc(jk))*c.rho_water / prhofirn(jk) ;
    % heat capacity in J/K
    heat_capa = heat_capa_vol .* volume;
    % energy needed to warm layer to 0degC 
    cold_content = deltaT*heat_capa; % in J
    warming_energy = min(cold_content, zsnout*c.L_fus*c.rho_water); % in J
    % warming layer
    ptsoil(jk) = ptsoil(jk) + warming_energy/heat_capa;
    % removing the corresponding amount from the prescribed melt
    zsnout = zsnout - warming_energy/c.L_fus/c.rho_water; % in m weq

    % Now the melting layer is either at melting point or zsnout is
    % depleted. We can use the rest of zsnout to change phase.
    
    %  How much frozen mass do we have in layer available for melting?       
    zdel = min(zsnout, psnowc(jk) + psnic(jk));

    % Update BV2017: Snow and ice are now melted simultaneously
    snow_melt = psnowc(jk)./(psnowc(jk)+psnic(jk)) .* zdel;
    ice_melt = psnic(jk)./(psnowc(jk)+psnic(jk)) .* zdel;
    
    psnowc(jk) = psnowc(jk) - snow_melt;
    psnic(jk) = psnic(jk) - ice_melt;
    pslwc(jk) = pslwc(jk) + zdel;
    
    zsnout = zsnout - zdel;
end

if sum((psnowc+psnic)==0)>1
    disp('MELTING MORE THAN ONE LAYER')
    disp(zsnmel)
end

end
