function   [QE,QH, T2m, q_2m, ws_10m, Ri] ...
            = SensLatFluxes_profile (z1,z2,t1,t2,q1,q2,u1,u2,p, c)
        
    % The virtual temperature is the temperature that dry dry air 
    % would have if its pressure and density were equal to those 
    % of a given sample of moist air (http://amsglossary.allenpress.com/glossary/)
    Tv1=  t1*(1.+0.61*q1);
    Tv2=  t2*(1.+0.61*q1);
    Tv_avg=(Tv1 + Tv2)/2.;

    % The potential temperature  is the temperature that an unsaturated
    % parcel of dry air would have if brought adiabatically and 
    % reversibly from its initial state to a standard pressure, p0,
    % typically 1000 hPa
    % Here we calculate the virtual potential temperature (dry)
    p0 = 1000; %standard pressure for the calculation of potential temperature

    theta_v1=t_v1*(p0/p)^c.kappa_poisson;
    theta_v2=t_v2*(p0/p)^c.kappa_poisson;
    dtheta = theta_v2-theta_v1;
    du = u2-u1;
    dz = z2-z1;
    z_r = dz/log(z2/z1);
    % z_r = sqrt(z1*z2);

    % Richardson number
    Ri = c.g ./ Tv_avg .* dtheta ./ dz .* (du./dz)^(-2);

    Ri(dudz <= 0 ) = NaN;
    if abs(Ri) > 2
        Ri = NaN;
    end
    
    % Stability criteria
    stab = stability(Ri);

    % air density
    ro_air =  p *100 /c.R_d/ Tv_avg ;

    [QH] = sensible_heat_flux (z_r,stab, Ri, ro_air,c);
    [QE]=  latent_heat_flux(z1,z2,q1,q2,u1,u2, stab, Ri, ro_air, c);

    QH(Ri >= 0.2 && Ri < 999) = 0;
    QE(Ri >= 0.2 && Ri < 999) = 0;
    
    T2m = t2;
    q_2m = q2;
    ws_10m = u2;
    end

    
%% Original IDL script
% pro Richardson,z1,z2,t1,t2,q2,q1,u1,u2,Ri
% 
% g=9.80665
% kk=273.15
% p=980.
% p0=1000.
% 
% Tv1=(T1+kk)*(1.+0.61*q1/1000.)
% Tv2=(T2+kk)*(1.+0.61*q2/1000.)
% 
% Tv_avg=(Tv1+Tv2)/2.
% 
% VPT1=Tv1*(p0/p)^0.286
% VPT2=Tv2*(p0/p)^0.286
% 
% DVPT=VPT2-VPT1
% 
% du=u2-u1
% 
% dz=z2-z1
% 
% dudz=du/dz
% 
% Ri=g/Tv_avg*DVPT/dudz^2.
% 
% invalid=where(u1 eq 999 or u2 eq 999,c) & if c gt 0 then Ri(invalid)=999.
% invalid=where(dudz le 0,c) & if c gt 0 then Ri(invalid)=999.
% invalid=where(dudz eq 999,c) & if c gt 0 then Ri(invalid)=999.
% invalid=where(dudz le 0,c)  & if c gt 0 then Ri(invalid)=999.
% invalid=where(abs(Ri) gt 2,c) & if c gt 0 then Ri(invalid)=999.
% 
% end
