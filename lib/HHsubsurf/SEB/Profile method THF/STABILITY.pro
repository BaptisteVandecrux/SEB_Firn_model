pro stability,Ri,stab

Ri_temp=Ri

b1=5.
b2=16.
a1=2.
a2=0.75

CRITICAL=WHERE(Ri_temp GE 0.2 and Ri_temp lt 999,C)
IF C GT 0 THEN Ri_temp(CRITICAL)=0.2

TOO_UNSTABLE=WHERE(Ri_temp LE -0.1,C)
IF C GT 0 THEN Ri_temp(TOO_UNSTABLE)=-1.

STAB=FLTARR(N_ELEMENTS(Ri_temp)) & STAB(*)=999.
; -------------------------------- STABLE CASE
STABLE=WHERE(Ri_temp GE 0. AND Ri_temp LE 0.2,C)
stab(STABLE)=(1-b1*Ri_temp(STABLE))^a1
; -------------------------------- UNSTABLE CASE
UNSTABLE=where(Ri_temp Ge -1. and Ri_temp LT 0.)
stab(UNSTABLE)=(1-b2*Ri_temp(UNSTABLE))^a2

invalid=where(Ri_temp eq 999.,c)
IF C GT 0 THEN stab(invalid)=999.

end