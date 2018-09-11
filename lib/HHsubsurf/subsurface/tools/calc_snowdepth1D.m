function [psn] =  calc_snowdepth1D (psnowc, psnic, snowbkt, c)
% calc_snowdepth1D: Diagnose snow depth (for use in albefor 
% parameterizations and other places?). Include only the snow above first
% perched ice layer (and include the snow of this layer)
%
%   This script was originally developped by Peter Langen (pla@dmi.dk) and
%   Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
%   Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================
notice = 1; % Should be read as "this layer is not ice"

% Update BV2017: including snowbkt
if ( snowbkt < c.smallno )
    psn = 0;
    % here we allow to look at the rest of the layers
else
    psn = snowbkt;
end

jk = 1;
% Then the next layers
while and(notice == 1, jk <= c.jpgrnd)
    if (psnic(jk) > c.icemax*c.cdel)
        %then it is considered as ice layer
        notice = 0;
    else
        psn = psn + psnowc(jk);
    end
    jk = jk+1;
end
end
