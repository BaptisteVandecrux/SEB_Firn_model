function [rho_snow] = IniRhoSnow(T, WS, elev, c)
%IniRhoSnow: Initialization of fresh snow density using different
%parametreization. If modified, check that the chosen parametrization is
%sent properly to the subsurface scheme.
%
% Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================

% INITIALIZATION OF SNOW DENSITY 
rho_snow = zeros(size(T));
mean_T = zeros(c.elev_bins);
for j=1:c.elev_bins
    mean_T = mean(T(:,j));
end

for j=1:c.elev_bins
%     mean_T(j) = mean(T(j,:));
    switch c.rho_snow_scheme
        case 0
            rho_snow(:,j) = c.fresh_snow_dens; %constant in time and place
        case 1
            
            rho_snow(:,j) = 625. + 18.7*(mean_T(j)-c.T_0) + 0.293*(mean_T(j)-c.T_0)^2;
            % surface snow density by Reeh et al. 2005
        case 2
            rho_snow(:,j) = 481. + 4.834*(mean_T(j)-c.T_0);
            % surface snow density by Kuipers-Munneke et al. 2015
        case 3
            rho_snow(:,j) = 369.58 -0.02298*elev(j); 
            % PLA RFO fresh snow density (May 2015): BV2017 using
            % elevation-based parametrization
        case 4
            rho_snow(:,j) = 350.; %default value everywhere
            ind = (T(:,j) >= 258.16);
            rho_snow(ind,j) = 50 + 1.7 * (T(ind,j) - 258.16).^1.5;
            clear ind
            ind = and((T(:,j) >= 258.16) , (WS(:,j) >= 5));
            rho_snow(ind,j) = 50 + 1.7 * (T(ind,j) - 258.16).^1.5 +...
                25 + 250 *(1-exp(-0.2*(WS(ind) - 5)));
            % surface snow density by Liston et al., 2007
            % only works above -15 degC
    end
end
end
