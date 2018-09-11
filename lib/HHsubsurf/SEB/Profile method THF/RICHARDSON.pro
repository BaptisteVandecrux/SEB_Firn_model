pro Richardson,z1,z2,t1,t2,q2,q1,u1,u2,Ri

g=9.80665
kk=273.15
p=980.
p0=1000.

Tv1=(T1+kk)*(1.+0.61*q1/1000.)
Tv2=(T2+kk)*(1.+0.61*q2/1000.)

Tv_avg=(Tv1+Tv2)/2.

VPT1=Tv1*(p0/p)^0.286
VPT2=Tv2*(p0/p)^0.286

DVPT=VPT2-VPT1

du=u2-u1

dz=z2-z1

dudz=du/dz

Ri=g/Tv_avg*DVPT/dudz^2.

invalid=where(u1 eq 999 or u2 eq 999,c) & if c gt 0 then Ri(invalid)=999.
invalid=where(dudz le 0,c) & if c gt 0 then Ri(invalid)=999.
invalid=where(dudz eq 999,c) & if c gt 0 then Ri(invalid)=999.
invalid=where(dudz le 0,c)  & if c gt 0 then Ri(invalid)=999.
invalid=where(abs(Ri) gt 2,c) & if c gt 0 then Ri(invalid)=999.

end