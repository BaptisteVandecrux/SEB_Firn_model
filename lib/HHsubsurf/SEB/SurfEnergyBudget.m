function [meltflux, Tsurf, dTsurf, EB_prev, stop] ...
    = SurfEnergyBudget (SRnet, LRin, Tsurf, k_eff, thick_first_lay, T_ice, T_rain,...
    dTsurf, EB_prev, SHF, LHF, rainfall, c)
% SurfEnergyBudget: calculates the surface temperature (Tsurf) and meltflux
% from the different elements of the energy balance. The surface
% temperature is adjusted iteratively until equilibrium between fluxes is
% found.
%
% Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================
    stop =0;
    % SURFACE ENERGY BUDGET -----------------------------------------------------------------------

    meltflux = SRnet(1) - SRnet(2) + LRin ...
        - c.em * c.sigma * Tsurf.^4 - (1 - c.em) * LRin ...
        + SHF + LHF ...
        -(k_eff(1)) * (Tsurf- T_ice(2)) / thick_first_lay ...
        + c.rho_water * c.c_w(1) * rainfall * c.dev ...
        / c.dt_obs *( T_rain - c.T_0);
        
    if meltflux >= 0 && Tsurf == c.T_0
        % stop iteration for melting surface
        stop =1;     
        return
    end
    
    if abs(meltflux) < c.EB_max
        % stop iteration when energy components in balance
        stop =1;   
        meltflux = 0;
        return
    end


    if meltflux/EB_prev < 0
        dTsurf = 0.5*dTsurf ;
        % make surface temperature step smaller when it overshoots EB=0
    end
    EB_prev = meltflux;

    %Update BV
    if meltflux < 0
        Tsurf = Tsurf - dTsurf ;
    else          
        Tsurf = min(c.T_0,Tsurf + dTsurf);
    end
end
