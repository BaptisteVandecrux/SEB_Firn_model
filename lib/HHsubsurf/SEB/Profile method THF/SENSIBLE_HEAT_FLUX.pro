; ------------------------------------------- Sensible Heat Flux
pro sensible_heat_flux,Ri,ro_air,Cp,k,z1,z2,t1,t2,q1,q2,u1,u2,QH,stab

moverr=3.4897
g=9.80665
kk=273.15
p=980.
p0=1000.

Richardson,z1,z2,t1,t2,q2,q1,u1,u2,Ri

stability,Ri,stab

; ----------------------------- du/dz
Tv1=(T1+kk)*(1.+0.61*q1/1000.)
Tv2=(T2+kk)*(1.+0.61*q2/1000.)

Tv_avg=(Tv1+Tv2)/2.

VPT1=Tv1*(p0/p)^0.286
VPT2=Tv2*(p0/p)^0.286

DVPT=VPT2-VPT1

du=u2-u1
dz=z2-z1

dudz=du/dz

invalid=where(u1 eq 999 or u2 eq 999,c)
if c gt 0 then dudz(invalid)=999.
if c gt 0 then du(invalid)=999.

valid=where(Tv_avg ne 999.)
ro_air=FLTARR(N_ELEMENTS(Ri)) & ro_air(*)=999.
ro_air(valid)=((moverr*(980.*100.))/(Tv_avg(valid)))/1000.

; -------------------------- choose lower profile
z_geometric=sqrt(z1*z2)


D_Theta_V=DVPT

valid=where(du ne 999. and stab ne 999.,c)

QH=FLTARR(N_ELEMENTS(Ri)) & QH(*)=999.

if c gt 0 then begin
  QH(valid)=ro_air(valid)*Cp*k^2*z_geometric^2*D_Theta_V(valid)*du(valid)/dz^2
endif

critical=where(Ri ge 0.2 and Ri lt 999,c) & if c gt 0 then QH(critical)=0.

x=where(abs(qh) gt 450,c) & if c gt 0 then QH(x)=999.


end