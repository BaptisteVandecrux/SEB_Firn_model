function Ri = richardson(z1,z2,T1,T2,q1,q2,u1,u2,p,c)
% function giving the richardson number developped by Jason Box
% translated from IDL by Baptiste Vandecrux
% b.vandecrux@gmail.com

    p0 = 1000; %standard pressure for the calculation of potential temperature
    
    % The virtual temperature is the temperature that dry dry air 
    % would have if its pressure and density were equal to those 
    % of a given sample of moist air (http://amsglossary.allenpress.com/glossary/)
    Tv1 = T1 .* (1.+0.61*q1/1000.); % virtual temperature
    Tv2 = T2 .* (1.+0.61*q2/1000.);
    Tv_avg = (Tv1+Tv2)/2;

    % The potential temperature  is the temperature that an unsaturated
    % parcel of dry air would have if brought adiabatically and 
    % reversibly from its initial state to a standard pressure, p0,
    % typically 1000 hPa
    VPT1 = Tv1*(p0/p)^c.kappa_poisson; %virtual potential temperatures
    VPT2 = Tv2*(p0/p)^c.kappa_poisson;

    DVPT=VPT2-VPT1;
    du = u2-u1;
    dz = z2-z1;
    
    Ri = c.g / Tv_avg * DVPT / dz * (du/dz)^(-2);

    Ri(dudz <= 0 ) = NaN;
    if abs(Ri) > 2
        Ri = NaN;
    end
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