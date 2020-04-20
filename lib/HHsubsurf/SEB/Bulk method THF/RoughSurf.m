function [z_h, z_q, u_star, Re] = RoughSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)

u_star = c.kappa * WS / (log(z_WS/z_0) - psi_m2 + psi_m1)  ;

% rough surfaces: Smeets & Van den Broeke 2008
Re = u_star * z_0 / nu;       

z_h = z_0 * exp(1.5 - 0.2*log(Re) - 0.11*(log(Re))^2); 

if z_h < 1e-6
     z_h = 1e-6;
end
z_q = z_h;

end
