function [hHirF] = hHirF(Theta,d, c)

%     ! The F in the name hHirF stands for "function"
% 
%     ! Calcuc.lates hydraulic suction h (in m) ac.ccording to 
%     ! Hirashima et al 2010.

alpha = 7.3*exp(1.9); %     ! Hirashima (15)
n = nHirF(d);
m=1-1/n; %               ! Hirashima (9)

Theta_nozero = max(Theta,c.smallno); % ! To avoid divide-by-zero
hHirF = 1/alpha * (Theta_nozero^(-1/m)-1)^(1/n); % ! Hirashima (9)
end
