function RH_wrt_w = RHice2water(RH_wrt_i, T, pres)
% function converting humidity with regards to ice into humidity with
% regards to water using GOFF-GRATCH (1945) formula.
%   
%   RH_wrt_i is an single value or an array of relative humidity with
%   regards to ice given in percent
%
%   T is the corresponding air temperature in degree Celsius
%
%   pres is the corresponding air pressure in hPa (not used in the current
%   form of the function)


T_0 = 273.15;
T_100 = 373.15;
% eps is the ratio of the molecular weights of water and dry air
eps = 0.62198; 

% Rv = 461.5;
% es_wtr = eps * exp( 1/Rv * ((L + T_0 * beta)*(1/T_0 - 1/T) - beta* log(T./T_0)));

% GOFF-GRATCH 1945 equation
es_wtr = 10.^( ...
    -7.90298*(T_100./T - 1) + 5.02808 * log10(T_100./T) ...   % saturation vapour pressure above 0 C (hPa)
    - 1.3816E-7 * (10.^(11.344*(1.-T./T_100))-1.) ...
    + 8.1328E-3*(10.^(-3.49149*(T_100./T-1)) -1.) + log10(1013.246) );

% WMO 2012 equation (based on Goff 1957)
% es_wtr = 10.^( ...
%     10.79574*(1 - T_0./T) + 5.028 * log10(T / T_0) ...   % saturation vapour pressure above 0 C (hPa)
%     + 1.50475E-4 * (1 - 10.^(-8.2969 * (T./T_0 - 1))) ...
%     + 0.42873E-3*(10.^(4.76955*(1 - T_0./T)) -1.) +  0.78614 + 2.0 );

es_ice = 10.^( ...
    -9.09718 * (T_0 ./ T - 1.) - 3.56654 * log10(T_0 ./ T) + ...
    0.876793 * (1. - T ./ T_0) + log10(6.1071)  );   % saturation vapour pressure below 0 C (hPa)

% es_ice = 10.^( ...
%     -9.09685 * (T_0 ./ T - 1.) - 3.56654 * log10(T_0 ./ T) + ...
%     0.87682 * (1. - T ./ T_0) + 0.78614 + 2.0  ); 

% q_sat_wtr = eps * es_wtr./(pres-(1-eps)*es_wtr);    % specific humidity at saturation over water
% q_sat = eps * es_ice./(pres-(1-eps)*es_ice);        % specific humidity at saturation over ice
freezing=ones(size(T));
freezing(T>=T_0) = 0;
freezing = freezing==1;

RH_wrt_w = RH_wrt_i;
RH_wrt_w(freezing) = RH_wrt_i(freezing) .*es_ice(freezing) ./ es_wtr (freezing) ;              % specific humidity in kg/kg
RH_wrt_w(RH_wrt_w<0)=0;
RH_wrt_w(RH_wrt_w>100) = 100;

end
