function [QE] = latent_heat_flux(z_r, phi_m, Ri, ro_air, c)
% function calculating the Latent Heat Flux QE as developped by Jason Box
% translated from IDL by Baptiste Vandecrux
% b.vandecrux@gmail.com
% ------------------------------------------- 
    QE = ro_air * c.L_vap * c.kappa^2 * z_r^2 * dq/dz * du/dz / phi_m.^2;
end

%% Original IDL script
% % ------------------------------------------- Sensible Heat Flux
% pro latent_heat_flux,Ri,ro_air,Cp,k,z1,z2,t1,t2,q1,q2,u1,u2,QE,stab
% 
% moverr=3.4897
% g=9.80665
% kk=273.15
% p=980.
% p0=1000.
% Lv=2.5e6
% 
% Richardson,z1,z2,t1,t2,q2,q1,u1,u2,Ri
% 
% stability,Ri,stab
% 
% % ----------------------------- du/dz
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
% dz=z2-z1
% 
% dudz=du/dz
% 
% valid=where(q1 ~= 999 and q2 ~= 999,c)
% dq=FLTARR(N_ELEMENTS(Ri)) & dq(*)=999.
% if c > 0 then dq(valid)=(q2(valid)-q1(valid))/1000.
% 
% invalid=where(u1 eq 999 or u2 eq 999,c)
% if c > 0 then dudz(invalid)=999.
% if c > 0 then du(invalid)=999.
% 
% valid=where(Tv_avg ~= 999.)
% ro_air=FLTARR(N_ELEMENTS(Ri)) & ro_air(*)=999.
% ro_air(valid)=((moverr*(980.*100.))/(Tv_avg(valid)))/1000.
% 
% % -------------------------- choose lower profile
% z_geometric=sqrt(z1*z2)
% 
% 
% D_Theta_V=DVPT
% 
% valid=where(du ~= 999. and stab ~= 999. and dq ~= 999,c)
% 
% QE=FLTARR(N_ELEMENTS(Ri)) & QE(*)=999.
% 
% if c > 0 then begin
%   QE(valid)=ro_air(valid)*Lv*k^2*z_geometric^2*dq(valid)*du(valid)/dz^2
% endif
% 
% critical=where(Ri >= 0.2 and Ri < 999,c) & if c > 0 then QE(critical)=0.
% 
% x=where(abs(QE) > 450,c) & if c > 0 then QE(x)=999.
% 
% 
% end