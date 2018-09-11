function [elev, pres, T, RH, WS, SRin, LRin, rho_atm, nu, mu, sum_dH_surf, Tdeep] ...
    = PrepareTransect(pres, T, RH, WS, SRin, LRin,c)
%PrepareTransect: Extrapolates the AWS data to each location of a specified
%transect.
%
% Author: Dirk Van As (dva@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================

% EXTRAPOLATING AWS DATA OVER THE TRANSECT 
if c.elev_bins ==  1
    c.elev_start = c.elev_AWS;
end
elev = (0:c.elev_bins-1)*c.elev_binsize + c.elev_start;

pres(1:c.M,:) = pres(1:c.M,1).*(T(1:c.M,1)./...
    (T(1:c.M,1)+(c.gradT+c.gradT_Tdep.*(T(1:c.M,1)-c.T_0)).*(elev-c.elev_AWS)))...
    .^(c.g./c.R_d./(c.gradT+c.gradT_Tdep.*(T(1:c.M,1)-c.T_0)));
T(1:c.M,:)  = T(1:c.M,1) + (elev-c.elev_AWS).*...
    (c.gradT + c.gradT_Tdep*(T(1:c.M,1)-c.T_0));

Tdeep = zeros(length(c.elev_bins));
if isempty(c.Tdeep_AWS)
    disp('Deep firn temperature not available, using the mean air temp instead.')
    c.Tdeep_AWS = nanmean(T);
end
for j = 1:c.elev_bins
    Tdeep(j) = c.T_0 + c.Tdeep_AWS; % longterm mean air temperature for KANU?
end

RH(1:c.M,:) = RH(1:c.M,1) + (elev-c.elev_AWS).*c.gradRH;
WS(1:c.M,:) = WS(1:c.M,1) + (elev-c.elev_AWS).*c.gradWS;
SRin(1:c.M,:) = SRin(1:c.M,1)*exp(-c.ext_air*(c.elev_AWS-elev));   % this needs inclusion of cloud effect and air mass dependence
LRin(1:c.M,:) = LRin(1:c.M,1) + (elev-c.elev_AWS)*(c.gradLRin + c.gradLRin_Tdep*(T(1:c.M,1)-c.T_0));

RHtoolow = find(RH < c.RH_min);
if sum(RHtoolow) > 0
    RH(RHtoolow) = c.RH_min;   % setting a lower limit for relative humidity (for supersaturation: see below)
end
WStoolow = find(WS < 0);
if sum(WStoolow) > 0
    WS(WStoolow) = 0  ;    % setting a lower limit for wind speed
end
rho_atm = 100*pres./c.R_d./T ;               % atmospheric density
mu=zeros(c.M,c.elev_bins);
nu=zeros(c.M,c.elev_bins);
mu = 18.27e-6.*(291.15+120)./(T+120).*(T./291.15).^1.5 ;  % dynamic viscosity of air (Pa s) (Sutherlands' formula using C = 120 K)
nu = mu./rho_atm  ;
% kinematic viscosity of air (m^2/s)
sum_dH_surf=0;

end
