  % ###########################################
  % PLA Darcy (Aug 2016)
  function [liqmax] =  CLliqF(rhosin, c)
  % ###########################################

  % Calculates per pore space volume irreducible liquid water content
  % following the parameterization by Coleou and Lesaffre (1998)

    % Cap input rho_snow at pore-close off value to avoid infinite numbers
    rhos = min(rhosin,c.rho_pco);
    % Porosity
    P  = 1-rhos/c.rho_ice;
    % Per mass parameterized Irreducible LWC
    Wm = 0.057 *  P/max(1-P, c.smallno) + 0.017;
    % Per pore space volume
    liqmax = Wm / ( P/max(1-P, c.smallno) ) ...
        * c.rho_ice/c.rho_water / max(1-Wm, c.smallno);
  end         