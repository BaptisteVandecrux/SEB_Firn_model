function [prhofirn, psnowc , psnic, pslwc, ptsoil , pdgrain] =...
    merge_small_layers (prhofirn, psnowc , psnic, pslwc, ptsoil , ...
        pdgrain, c)

% Now that the water has percolated and refrozen, we go through the column and check
% that there is no layer that ended up too small

for jk =1:c.jpgrnd-1

    too_small = (psnic(jk) + psnowc(jk) + pslwc(jk)) < c.lim_new_lay;
    while too_small== 1
        if (psnic(jk) + psnowc(jk) + pslwc(jk)) > c.smallno
            % then we merge this layer with the layer below
            if (psnowc(jk+1) + psnowc(jk))~=0;
                pdgrain(jk+1) = (psnowc(jk) * pdgrain(jk) + psnowc(jk+1) * pdgrain(jk+1)) ...
                   /(psnowc(jk+1) + psnowc(jk));
                prhofirn(jk+1) = (psnowc(jk+1) + psnowc(jk)) ...
                   / (psnowc(jk)/prhofirn(jk) +  psnowc(jk+1)/prhofirn(jk+1));
            else
                pdgrain(jk+1) = c.dgrainNew;
                prhofirn(jk+1) = c.rho_ice;
            end
            ptsoil(jk+1) = ((psnic(jk) + psnowc(jk) + pslwc(jk)) * ptsoil(jk) ...
                + (psnic(jk+1) + psnowc(jk+1) + pslwc(jk+1)) * ptsoil(jk+1)) ...
                /(psnic(jk) + psnowc(jk) + pslwc(jk) + psnic(jk+1) + psnowc(jk+1) + pslwc(jk+1));

            psnowc(jk+1) = psnowc(jk+1) + psnowc(jk);
            psnic(jk+1) = psnic(jk+1) + psnic(jk);
            pslwc(jk+1) = pslwc(jk+1) + pslwc(jk);

            % now layer jk is free but we need to shift the column from 1 to jk
            % in order to free the first layer instead (just for being
            % compatible with the SplitLayer function)
            if jk>1
                psnowc(2:jk) = psnowc(1:jk-1);
                psnic(2:jk) = psnic(1:jk-1);
                pslwc(2:jk) = pslwc(1:jk-1);
                pdgrain(2:jk) = pdgrain(1:jk-1);
                prhofirn(2:jk) = prhofirn(1:jk-1);
                ptsoil(2:jk) = ptsoil(1:jk-1);
                % now top layer is free and another layer can be splitted
                % elsewhere
            end
        end

        [psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil] = ...
            split_layer(psnic, psnowc, pslwc, pdgrain, prhofirn, ptsoil, c);

        % if a layer was split above indice jk then the new, merged layer
        % is now located at jk+1 and will be the next one tested in the
        % for loop. However if the split happened below indice jk then
        % the merged layer is still located at indice jk and ti should be
        % tested again to see if this merging was enough
        too_small = (psnic(jk) + psnowc(jk) + pslwc(jk)) < c.lim_new_lay;
        
        % in the first case, the new 'too_small' will be the one of the layer
        % just above the one that was just merged and should be 0 because
        % the loop already went through that layer. In the second case, the
        % new 'too_small' is the one of the merged layer and if still too
        % small, would be set to 1 and another layer would be merged to it.
    end
end
end