function [slwc] = MimicAquiferFlow(snowc, rhofirn, snic, slwc, k,  c)
% here is a small option that allows to inject water in the
% subsurface to mimic the inflow from the aquifer
            % flow from aquifer 
            % Poinar et al. (2017): 10000 m3/yr = 11.4 m weq / hr
            % Poinar et al. (2017): 5000 m3/yr = 0.6 m weq / hr
            % trial: 1000 m3/yr = 0.11 m weq / hr
            % trial: 500 m3/yr = 0.006 m weq / hr
            thickness_act = snowc.*(c.rho_water./rhofirn) + snic.*(c.rho_water./c.rho_ice);

            depth_act=zeros(size(thickness_act));
            for jj=1:length(thickness_act(1,:))
                    depth_act(1:end,jj)=cumsum(thickness_act(1:end,jj));
            end

            [~, ind_15] = min(abs(depth_act-15));
            [~, ind_25] = min(abs(depth_act-25));
            thick_aq = depth_act(ind_25)- depth_act(ind_15);
            
%             pore_space = snowc.*c.rho_water .*(c.rho_ice - rhofirn)./rhofirn./c.rho_ice;
%             pore_space(snowc<c.smallno) = 0;
            if k>1
                slwc(ind_15:ind_25) = slwc(ind_15:ind_25) + 11.4*thickness_act(ind_15:ind_25)/thick_aq;
            end
            % surface flow from catchment
            slwc(1:5) = slwc(1:5) + snmel*thickness_act(1:5)/(depth_act(5)- depth_act(1));
            
end