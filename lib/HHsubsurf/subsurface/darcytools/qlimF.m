function [qlimF] = qlimF(pl1,pl2,ps1,ps2,pi1,pi2,rhos1,rhos2,d1,d2, c)
% The F in the name stands for "function"
% Calculates Hirashima (2010) qlim in eqn 20  iteratively

% Physical layer thicknesses (m)
delz1 = c.rho_water*(ps1/rhos1 + pi1/c.rho_ice); % layer 1
delz2 = c.rho_water*(ps2/rhos2 + pi2/c.rho_ice); % layer 2
delz = (delz1+delz2)/2;             % Midpoint-to-midpoint distance

% First, see if there is a solution, by moving all in upper layer down.
% If this provides positive h1-(h2+delz) then there is a solution
qlim = pl1;

pl1test = pl1-qlim;
pl2test = pl2+qlim;
Theta1 = ThetaF(pl1test,ps1,rhos1, c);
Theta2 = ThetaF(pl2test,ps2,rhos2, c);
h1 = hHirF(Theta1,d1, c);
h2 = hHirF(Theta2,d2, c);
diff = h1 - (h2 + delz);
if ( diff < 0 )
    Conv = 1;  % If moving everything isn't enough for equilibrium, move everything and do no more
    qlimOut = pl1;
else
    Conv = 0;  % If it is enough, we haven't converged yet and we start iteration below
end

if ( Conv == 0 )
    % First guess is half of water in top layer
    qlim = pl1/2;
    qlimL =0;
    qlimR =pl1;
    
    Nmax = 11; % Number of iterations (usually 10 suffices since
    % we are halving intervals and 0.5^9 - 0.5^10 = 9.7656e-04 < 0.001)
    
    for i = 1:Nmax
        pl1test = pl1-qlim;
        pl2test = pl2+qlim;
        Theta1 = ThetaF(pl1test,ps1,rhos1, c);
        Theta2 = ThetaF(pl2test,ps2,rhos2, c);
        h1 = hHirF(Theta1,d1, c);
        h2 = hHirF(Theta2,d2, c);
        diff = h1 - (h2 + delz); % Difference between LHS and RHS in Hirashima eqn 20
        
        % If positive, we moved too much, and  qlim becomes new right end point
        % if negative, we moved too little, and qlim becoms new left end:
        if (diff > 0 )
            qlimR = qlim;
        else
            qlimL = qlim;
        end
        % New value is halfway between new interval end points
        qlim = (qlimR+qlimL)/2;
    end
    % Set final output qlim to iterated value:
    qlimOut = qlim;
end

qlimF = qlimOut;
end
