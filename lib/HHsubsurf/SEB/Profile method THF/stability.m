function phi_m = stability(Ri)
% function calculating the correction parameters for the similarity 
% developped by Jason Box
% translated from IDL by Baptiste Vandecrux
% b.vandecrux@gmail.com

    Ri_temp=Ri;

    critical = and(Ri_temp >= 0.2 , Ri_temp < 999);
    if sum(critical) > 0 
        Ri_temp(critical) = 0.2;
    end

    too_unstable = (Ri_temp <= -0.1);
    if sum(too_unstable) > 0
        Ri_temp(too_unstable) = -1;
    end

    % -------------------------------- STABLE CASE
    phi_m = NaN(size(Ri_temp));
    stable = and(Ri_temp >= 0 , Ri_temp <= 0.2);
    % according to Box and Steffen (2001)
    phi_m(stable) = (1-5.2*Ri_temp(stable))^(-0.5);  

    % -------------------------------- UNSTABLE CASE
    unstable = and(Ri_temp >= -1, Ri_temp < 0);
    phi_m(unstable) = (1-18*Ri_temp(unstable))^(-0.25); % Box and Steffen (2001)
end

%% Original IDL script
% pro stability,Ri,stab
% 
% Ri_temp=Ri
% 
% b1=5.
% b2=16.
% a1=2.
% a2=0.75
% 
% CRITICAL=WHERE(Ri_temp GE 0.2 and Ri_temp lt 999,C)
% IF C > 0 THEN Ri_temp(CRITICAL)=0.2
% 
% TOO_UNSTABLE=WHERE(Ri_temp <= -0.1,C)
% IF C > 0 THEN Ri_temp(TOO_UNSTABLE)=-1.
% 
% STAB=FLTARR(N_ELEMENTS(Ri_temp)) & STAB(*)=999.
% ; -------------------------------- STABLE CASE
% STABLE=WHERE(Ri_temp GE 0. AND Ri_temp <= 0.2,C)
% stab(STABLE)=(1-b1*Ri_temp(STABLE))^a1
% ; -------------------------------- UNSTABLE CASE
% UNSTABLE=where(Ri_temp Ge -1. and Ri_temp LT 0.)
% stab(UNSTABLE)=(1-b2*Ri_temp(UNSTABLE))^a2
% 
% invalid=where(Ri_temp eq 999.,c)
% IF C > 0 THEN stab(invalid)=999.
% 
% end


