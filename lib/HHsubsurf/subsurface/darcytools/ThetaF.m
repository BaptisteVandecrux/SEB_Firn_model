function [ThetaF] = ThetaF(pl,ps,rhos, c)
%     ! The F in the name stands for "function"
%     ! Calculates effective water saturation.
%     ! Force Theta to be between 0 and 1
% Calculate liqmax according to Coleou and Lesaffre?
% liqmaxloc is the value of liqmax used locally here
if ( c.calc_CLliq )
  liqmaxloc = CLliqF(rhos, c);
else
  liqmaxloc = c.liqmax;
end

if (ps > c.smallno)
    %        ! Force Theta to be between 0 and 1
    ThetaF = min(1, ...
        max ( 0, ...
        ( (pl / ps * rhos * c.rho_ice / c.rho_water / (c.rho_ice - rhos)...
        - liqmaxloc) / (1-liqmaxloc) ) ));
else
    %        ! If there is no snow, write out Theta=1 to avoid divide-by-zero
    ThetaF = 1;
end
end

