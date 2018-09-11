function [zsn_capaF] = zsn_capaF(RHOS,ptsoil)
% The F in the name zsn_capaF stands for "function"
zsn_capaF = cpiceF(ptsoil).*RHOS;  % snow vol. heat capacity   [J/m**3/K]
% Here the specific (per kg) heat is converted into volumetric heat (per
% m^3) capacity by multiplying by the density.

end
