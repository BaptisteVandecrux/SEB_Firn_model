function   [QE,QH, T2m, q_2m, ws_10m, Ri] ...
            = SensLatFluxes_profile (z1,z2,t1,t2,q1,q2,u1,u2,p, c)
        
Ri = richardson(z1,z2,t1,t2,q1,q2,u1,u2,p,c);
stab = stability(Ri);

% The virtual temperature is the temperature that dry dry air 
% would have if its pressure and density were equal to those 
% of a given sample of moist air (http://amsglossary.allenpress.com/glossary/)
Tv1= t1*(1.+0.61*q1);
Tv2= t2*(1.+0.61*q2);

Tv_avg=(Tv1+Tv2)/2.;
% air density
ro_air =  p *100 /c.R_d/ Tv_avg ;

[QH]=  sensible_heat_flux (z1,z2,t1,t2,q1,q2,u1,u2,p, stab, Ri, ro_air, c);
[QE]=  latent_heat_flux(z1,z2,q1,q2,u1,u2, stab, Ri, ro_air, c);

T2m = t2;
q_2m = q2;
ws_10m = u2;
end