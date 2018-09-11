    
% ORIGINAL FORTRAN CODE
%      module darcytools
% 
%   implicit none
% 
% contains
%   ! ###########################################
%   ! PLA Darcy (May 2016)
%   function hHirF(Theta,d)
%     ! ###########################################
%     use declare, only : smallno
%     IMPLICIT NONE
%     real(8), intent(in)  :: Theta, d
%     real(8)              :: alpha, n, m, hHirF, Theta_nozero
%     ! The F in the name hHirF stands for "function"
% 
%     ! Calculates hydraulic suction h (in m) according to 
%     ! Hirashima et al 2010.
% 
%     alpha = 7.3*d+1.9     ! Hirashima (15)
%     n = nHirF(d)
%     m=1-1/n               ! Hirashima (9)
% 
%     Theta_nozero = max(Theta,smallno) ! To avoid divide-by-zero
%     hHirF = 1/alpha * (Theta_nozero**(-1/m)-1)**(1/n) ! Hirashima (9)
%   end function hHirF
% 
%   ! ###########################################
%   ! PLA Darcy (May 2016)
%   function nHirF(d)
%     ! ###########################################
%     IMPLICIT NONE
%     real(8), intent(in)  :: d
%     real(8)              :: nHirF
%     ! Calculates n as in Hirashima et al 2010.
%     nHirF = 15.68*exp(-.46*d)+1 ! Hirashima (17)
%   end function nHirF
% 
%   ! ###########################################
%   ! PLA Darcy (May 2016)
%   function ThetaF(pl,ps,rhos)
%     ! ###########################################
%     use declare, only : rhoh2o, rhoice, liqmax, smallno
%     IMPLICIT NONE
%     real(8), intent(in)  :: pl,ps,rhos
%     real(8)              :: ThetaF
% 
%     ! The F in the name stands for "function"
%     ! Calculates effective water saturation.
%     ! Force Theta to be between 0 and 1
%     if (ps .gt. smallno) then
%        ! Force Theta to be between 0 and 1
%        ThetaF = min(1.d0, max ( ( (pl/ps* rhos*rhoice/rhoh2o/(rhoice-rhos) - liqmax) / (1-liqmax) ) , 0.d0 ))
%     else
%        ! If there is no snow, write out Theta=1 to avoid divide-by-zero
%        ThetaF = 1.d0
%     end if
%   end function ThetaF
% 
%   ! ###########################################
%   ! PLA Darcy (May 2016)
%   function kF(Theta,d,rhos,pi,ps)
%   ! ###########################################
%     use declare, only : nuw, g_grav, rhoice, kice, whwice
%     IMPLICIT NONE
%     real(8), intent(in)  :: Theta, d, rhos, pi, ps
%     real(8)              :: kr, ks, k11, k22factor, Hs, Hi, fsnow, n, m
%     real(8)              :: kF
% 
%     ! The F in the name stands for "function"
%     ! Permeability as in Hirashima 2010, using units on Ks as in
%     ! Calonne (2012) and correction for ice layering using Colbeck 1975 (eqn 32)
% 
%     ! Saturated hydraulic conductivity of snow (m/s) as fitted by Shimizu (1970). Here
%     ! units are as in Calonne et al (2012) and d/1000 is the grain diameter in m. It
%     ! corresponds to H10-eqn3 with units as in Calonne
%     ks = g_grav/nuw*0.077*(d/1000.)**2. * exp(-7.8d-3*rhos) 
%     
%     ! Unsaturated correction to conductivity. Hirashima eqn 10
%     n = nHirF(d)
%     m=1-1/n               ! Hirashima (9)
%     kr = Theta**0.5 * ( 1.-(1.-Theta**(1./m))**m )**2. ! Hirashima eqn 10
% 
%     ! Unsaturated hydraulic conductivity of the snow. Hirashima eqn 11
%     k11 = kr*ks
% 
%     ! Factor to divide k11 by to account for ice lenses within the
%     ! layer. Colbeck 1975 eqn 32
%     Hs    = ps/rhos      ! total depth of snow in layer
%     Hi    = pi/rhoice    ! total depth of ice in layer
%     fsnow = Hs / (Hs+Hi) ! Fraction of layer that is snow
%     if (k11 .gt. 0.d0) then
%       k22factor = fsnow + (1.-fsnow) * (1.+whwice)/(kice/k11 + whwice)
%     else
%       ! If k11 is 0 (because Theta is 0) then set to 1, so kF = 0/1 = 0 in next line
%       k22factor = 1.d0
%     end if
%     ! Effective hydraulic conductivity (in vertical direction perpendicular to
%     ! ice lenses) in snow-ice-sublayered model. Colbeck 1975 eqn 32
%     kF = k11/k22factor
%   end function kF
% 
%   ! ###########################################
%   ! PLA Darcy (May 2016)
%   function qlimF(pl1,pl2,ps1,ps2,pi1,pi2,rhos1,rhos2,d1,d2)
%     ! ###########################################
%     use declare, only : rhoice, rhoh2o
%     IMPLICIT NONE
%     real(8), intent(in)  :: pl1, pl2, ps1, ps2, pi1, pi2, rhos1, rhos2, d1, d2
%     real(8)              :: delz, delz1, delz2, qlim, pl1test, pl2test, Theta1, Theta2
%     real(8)              :: h1, h2, diff, qlimL, qlimR, qlimOut
%     integer              :: i, Conv, Nmax
%     real(8)              :: qlimF
% 
%     ! The F in the name stands for "function"
%     ! Calculates Hirashima (2010) qlim in eqn 20  iteratively
% 
%     ! Physical layer thicknesses (m)
%     delz1 = rhoh2o*(ps1/rhos1 + pi1/rhoice) ! layer 1
%     delz2 = rhoh2o*(ps2/rhos2 + pi2/rhoice) ! layer 2
%     delz = (delz1+delz2)/2.             ! Midpoint-to-midpoint distance
%     
%     ! First, see if there is a solution, by moving all in upper layer down.
%     ! If this provides positive h1-(h2+delz) then there is a solution
%     qlim = pl1
% 
%     pl1test = pl1-qlim
%     pl2test = pl2+qlim
%     Theta1 = ThetaF(pl1test,ps1,rhos1)
%     Theta2 = ThetaF(pl2test,ps2,rhos2)
%     h1 = hHirF(Theta1,d1)
%     h2 = hHirF(Theta2,d2)
%     diff = h1 - (h2 + delz)
%     if ( diff .lt. 0. ) then
%        Conv = 1  ! If moving everything isn't enough for equilibrium, move everything and do no more
%        qlimOut = pl1
%     else
%        Conv = 0  ! If it is enough, we haven't converged yet and we start iteration below
%     endif
% 
%     if ( Conv .eq. 0 ) then
%        ! First guess is half of water in top layer    
%        qlim = pl1/2.
%        qlimL =0.
%        qlimR =pl1
% 
%        Nmax = 11 ! Number of iterations (usually 10 suffices since
%        ! we are halving intervals and 0.5^9 - 0.5^10 = 9.7656e-04 < 0.001)
% 
%        do i = 1,Nmax
%           pl1test = pl1-qlim
%           pl2test = pl2+qlim
%           Theta1 = ThetaF(pl1test,ps1,rhos1)
%           Theta2 = ThetaF(pl2test,ps2,rhos2)
%           h1 = hHirF(Theta1,d1)
%           h2 = hHirF(Theta2,d2)
%           diff = h1 - (h2 + delz) ! Difference between LHS and RHS in Hirashima eqn 20
% 
%           ! If positive, we moved too much, and  qlim becomes new right end point
%           ! if negative, we moved too little, and qlim becoms new left end:
%           if (diff .gt. 0 ) then
%              qlimR = qlim
%           else
%              qlimL = qlim
%           end if
%           ! New value is halfway between new interval end points
%           qlim = (qlimR+qlimL)/2.
%        end do
%        ! Set final output qlim to iterated value:
%        qlimOut = qlim
%     end if
% 
%     qlimF = qlimOut
%   end function qlimF
% 
%   ! ###########################################
%   ! PLA Darcy (May 2016)
%   function dgraindtF(dg,pl,ps)
%   ! ###########################################
%     use declare, only : pi_val, smallno
%     IMPLICIT NONE
%     real(8), intent(in)  :: dg, pl, ps
%     real(8)              :: L, Tusima, Brun
%     real(8)              :: dgraindtF 
% 
%     ! The F in the name stands for "function"
%     ! Calculates d(grainsize-diameter) /dt (in mm/s) as in Hirashima (2010), using saturated
%     ! and unsaturated rates from Brun (1989), Tusima (1978) and Katsuhima (2009).
% 
%     ! Mass liquid water content in %
%     L = pl/(ps+smallno)*100.d0
% 
%     ! Brun derives a function for grain growth that increases as L^3 (and divides by d^2).
%     ! Katushima uses this for L<=10%. Beyond 10 % Katushima uses the minimum value of the
%     ! Brun formulation and the constant-divided-by-d^2 found by Tusima.
%     ! However, the Tusima-constant is lower than the L^3 function at 10% so that would lead
%     ! to a discontinuous drop when going beyond 10%.
%     ! Instead, we do as we also think Hirashima does "Therefore, to consider grain growth
%     ! in the water-saturated layer, the grain coarsening model of Tusima (1978) was used
%     ! as an upper boundary to saturated grain growth following Katsushima et al. (2009a)":
% 
%     ! Take the L^3 function by Brun up until this becomes larger than the Tusima-constant:
%     Brun   = 2./pi_val * ( 1.28d-8 + 4.22d-10 * L**3.)
%     Tusima = 6.94d-8
%     ! d (grainsize-diamter) /dt in mm/s
%     dgraindtF = 1./(dg**2.)*min(Brun,Tusima) 
%   end function dgraindtF
% 
% end module darcytools

     