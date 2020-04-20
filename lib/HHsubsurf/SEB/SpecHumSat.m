function [RH, q] = SpecHumSat(RH, T, pres, c)

% SpecHumSat
% - calculates saturation vapour pressure of water (es_wtr) ice (es_ice)
% given the temperature T
% - corrects RH for supersaturation
% - calculates the specific humidity at saturation (qsat) for sub- and supra-
% freezing conditions
% - calculates the specific humidity as function of RH and qsat
%
% Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================


% SPECIFIC HUMIDITY & SATURATION -----------------------------------------------------------------------
es_wtr = 10.^(-7.90298*(c.T_100./T-1) ...
    + 5.02808 * log10(c.T_100./T) ...   % saturation vapour pressure above 0 C (hPa)
    - 1.3816E-7 * (10.^(11.344*(1.-T./c.T_100))-1.) ...
    + 8.1328E-3*(10.^(-3.49149*(c.T_100./T-1)) -1.) + log10(c.es_100));

es_ice = 10.^(-9.09718 * (c.T_0 ./ T - 1.) ...
    - 3.56654 * log10(c.T_0 ./ T) + ...
    0.876793 * (1. - T ./ c.T_0) + log10(c.es_0));   % saturation vapour pressure below 0 C (hPa)
q_sat = c.es * es_wtr./(pres-(1-c.es)*es_wtr);    % specific humidity at saturation (incorrect below melting point)

freezing = find(T < c.T_0);            % replacing saturation specific humidity values below melting point
if sum(freezing) > 0
    q_sat(freezing) = c.es * es_ice(freezing)./(pres(freezing)-(1-c.es)*es_ice(freezing));
end
% supersaturated = find(RH > 100) ;      % replacing values of supersaturation by saturation
% if sum(supersaturated) > 0
%     RH(supersaturated) = 100;
% end
q = RH.*q_sat/100  ;              % specific humidity in kg/kg
end

