function  [z_h, z_q, u_star, Re] = SmoothSurf(WS, z_0, psi_m1, psi_m2, nu, z_WS, c)

u_star = c.kappa*WS/(log(z_WS/z_0)- psi_m2 + psi_m1)  ;
Re = u_star * z_0 / nu; 

if Re <= 0.135
    range = 1;
elseif and(Re > 0.135 , Re < 2.5)
    range = 2;
elseif Re >= 2.5
    range = 3;
else
    disp('ERROR')
    disp(Re)
    disp(u_star)
    disp(z_WS)
    disp(psi_m2)
    disp(psi_m1)
    disp(nu)
    disp(range)
end

 % smooth surfaces: Andreas 1987
z_h =  z_0 * exp(c.ch1(range) + c.ch2(range)*log(Re) + c.ch3(range)*(log(Re))^2);  
z_q = z_0 * exp(c.cq1(range) + c.cq2(range)*log(Re) + c.cq3(range)*(log(Re))^2);

    if z_h < 1e-6
         z_h = 1e-6;
    end
    if z_q < 1e-6
         z_q = 1e-6;
    end
end
