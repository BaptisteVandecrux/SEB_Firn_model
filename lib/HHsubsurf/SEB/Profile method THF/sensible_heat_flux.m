function [QH] = sensible_heat_flux (z_r,phi_m, Ri, ro_air,c)
% function calculating the Sensible Heat Flux QH usingthe profile method
% as developped by Jason Box.
% translated from IDL by Baptiste Vandecrux
% b.vandecrux@gmail.com
% ------------------------------------------- 

% phi_h needs to be tuned up for very unstable conditions according to Box and Steffen (2001)
phi_m(Ri<-0.03) = phi_m(Ri<-0.03)*1.3;
QH =  ro_air .* c.c_pd .* c.kappa^2 .* z_r.^2 .* dtheta ./dz .* du./dz ./ phi_m.^2;
end

%% Original IDL script
% ; ------------------------------------------- Sensible Heat Flux
% pro sensible_heat_flux,Ri,ro_air,Cp,k,z1,z2,t1,t2,q1,q2,u1,u2,QH,phi_m
% 
% moverr=3.4897
% g=9.80665
% c.T_0=273.15
% p=980.
% p0=1000.
% 
% Richardson,z1,z2,t1,t2,q2,q1,u1,u2,Ri
% 
% phi_mility,Ri,phi_m
% 
% ; ----------------------------- du/dz
% Tv1=(T1+c.T_0)*(1.+0.61*q1/1000.)
% Tv2=(T2+c.T_0)*(1.+0.61*q2/1000.)
% 
% Tv_avg=(Tv1+Tv2)/2.
% 
% VPT1=Tv1*(p0/p)^0.286
% VPT2=Tv2*(p0/p)^0.286
% 
% DVPT=VPT2-VPT1
% 
% du=u2-u1
% dz=z2-z1
% 
% dudz=du/dz
% 
% invalid=where(u1 == 999 or u2 == 999,c)
% if c gt 0 then dudz(invalid)=999.
% if c gt 0 then du(invalid)=999.
% 
% valid=where(Tv_avg ne 999.)
% ro_air=FLTARR(N_ELEMENTS(Ri)) & ro_air(*)=999.
% ro_air(valid)=((moverr*(980.*100.))/(Tv_avg(valid)))/1000.
% 
% ; -------------------------- choose lower profile
% z_geometric=sqrt(z1*z2)
% 
% 
% D_Theta_V=DVPT
% 
% valid=where(du ne 999. and phi_m ne 999.,c)
% 
% QH=FLTARR(N_ELEMENTS(Ri)) & QH(*)=999.
% 
% if c gt 0 then begin
%   QH(valid)=ro_air(valid)*Cp*k^2*z_geometric^2*D_Theta_V(valid)*du(valid)/dz^2
% endif
% 
% critical=where(Ri ge 0.2 and Ri < 999,c) & if c gt 0 then QH(critical)=0.
% 
% x=where(abs(qh) gt 450,c) & if c gt 0 then QH(x)=999.
% 
% 
% end