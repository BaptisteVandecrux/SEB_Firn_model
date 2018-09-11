function [kF] = kF(Theta,d,rhos,pi,ps, c)
%     ! The F in the name stands for "function"
%     ! Permeability as in Hirashima 2010, using units on Ks as in
%     ! Calonne (2012) and correction for ice layering using Colbeck 1975 (eqn 32)

%     ! Saturated hydraulic conductivity of snow (m/s) as fitted by Shimizu (1970). Here
%     ! units are as in Calonne et al (2012) and d/1000 is the grain diameter in m. It
%     ! corresponds to H10-eqn3 with units as in Calonne
% ks_shim = c.g/c.nuw*0.077*(d/1000.)^2. * exp(-7.8e-3*rhos) ;
% Permeability according to Calonne (2012)
ks = 3*(d/2000)^2*c.g/c.nuw*exp(-0.013*rhos);

% fprintf('Permeability %f %f\n',ks/c.g*c.nuw/(d/2000)^2, ks_shim/c.g*c.nuw/(d/2000)^2);

% Unsaturated correction to conductivity. Hirashima eqn 10
n = nHirF(d);
m=1-1/n ;              % Hirashima (9)
kr = Theta^0.5 * ( 1.-(1.-Theta^(1./m))^m )^2; % Hirashima eqn 10

% Unsaturated hydraulic conductivity of the snow. Hirashima eqn 11
k11 = kr*ks;

% Factor to divide k11 by to account for ice lenses within the
% layer. Colbeck 1975 eqn 32
Hs    = ps/rhos ;     % total depth of snow in layer
Hi    = pi/c.rho_ice ;   % total depth of ice in layer
fsnow = Hs / (Hs+Hi); % Fraction of layer that is snow
if (k11 > 0)
    k22factor = fsnow + (1-fsnow) * (1.+c.whwice)/(c.kice/k11 + c.whwice);
else
    % If k11 is 0 (because Theta is 0) then set to 1, so kF = 0/1 = 0 in next line
    k22factor = 1;
end
% Effective hydraulic conductivity (in vertical direction perpendicular to
% ice lenses) in snow-ice-sublayered model. Colbeck 1975 eqn 32
kF = k11/k22factor;

end
