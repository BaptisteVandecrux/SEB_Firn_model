function [pgrndc,pgrndd,   pgrndcapc, pgrndhflx ]...
    = update_tempdiff_params (prhofirn, pTdeep                    ...
    , psnowc, psnic, ptsoil, zso_cond, zso_capa, c)
% update_tempdiff_params: Updates the thermal capacity and conductivity of 
% subsurface column based on the new density and temperature profiles. Also
% calcualates subsurface heat flux to the surface pgrndhflx and calorific capacity of the
% ground pgrndcapc.
% Input:
%   zso_cond - Thermal conductivity of pure ice. Calculated in ice_heats.
%   zso_capa - Specific heat capacity of ice and snow.
%   pgrndc - heat capacity (W/K)
% 
% This script was originally developped by Peter Langen (pla@dmi.dk) and
% Robert S. Fausto (rsf@geus.dk) in FORTRAN then translated to Matlab by 
% Baptiste Vandecrux (bava@byg.dtu.dk).
%=========================================================================


% PETER AND RUTH's VERSION (Note: we ignore the liquid content in all of 
% the following). We need physical layer thicknesses (i.e., in
% snow and ice meters, not liquid meters):

% We consider that each of our temperature in the column are taken at
% the center of the cell. In other words they are seen as nodes and not 
% as cells anymore. We then calculate the heat capacity and thermal
% conductivity of the material that separates neighbouring nodes.

snowV1=zeros(c.jpgrnd,1);
snowV2=zeros(c.jpgrnd,1);
iceV=zeros(c.jpgrnd,1);
pgrndc=zeros(c.jpgrnd,1);
pgrndd=zeros(c.jpgrnd,1);


%% Calculate midpoint volume-weighted versions of capa (heat capacity) and
% kappa (thermal conductivity)

% PLA densification (Feb 2015) Update BV2017
% The first layer should be in equilibrium with the surface temperature.
% Therefore first two nodes are seperated by all of layer 1 and half of
% layer 2.
% volume of snow in upper layer:
snowV1(1) =  psnowc(1)*c.rho_water/prhofirn(1); 
% volume of snow in lower half layer:
snowV2(1) = 0.5*psnowc(2)*c.rho_water/prhofirn(2); 
% volume of ice in both layer:
iceV(1)  = (psnic(1) + 0.5*psnic(2)) *c.rho_water / c.rho_ice; 

% For following nodes, two neighboring nodes are separated by half of the 
% upper layer and half of the lower layer
% volume of snow in upper half layer:
snowV1(2:c.jpgrnd) = 0.5*psnowc(2:c.jpgrnd)  .*...
    c.rho_water./prhofirn(2:c.jpgrnd);
% volume of snow in lower half layer:
snowV2(2:c.jpgrnd - 1) = 0.5*psnowc(3:c.jpgrnd)*c.rho_water./prhofirn(3:c.jpgrnd);
% volume of ice in both half layers:
iceV(2:c.jpgrnd - 1)  = 0.5*(psnic(2:c.jpgrnd - 1) + psnic(3:c.jpgrnd ))...
    .*c.rh2oice;

% Bottom layer zcapa assuming below is ice to same thickness (at least)
snowV2(c.jpgrnd) = 0;
iceV(c.jpgrnd) = 0.5 * (psnic(c.jpgrnd) + c.cdel(c.jpgrnd)) *c.rh2oice;

% total mass separating two nodes
totalV = snowV1 + snowV2 + iceV;

% ice and snow volumetric fractions in the material separating two nodes
snow_frac_lay_1 = snowV1./ totalV;
snow_frac_lay_2 = snowV2./ totalV;
ice_frac = iceV ./ totalV;

% heat capacity of the layer calculated as the volume-weighted average of
% the snow and ice heat capacity. Mind the repeated index for the last
% layer. Here zcapa is still volumetric since everything on the right hand
% side has been deivided by 'totalV'. It is in J/m^3/K.
zcapa = snow_frac_lay_1 .* zsn_capaF(prhofirn, ptsoil) +...
    snow_frac_lay_2 .* ...
    zsn_capaF(prhofirn([2:c.jpgrnd, c.jpgrnd]), ptsoil([2:c.jpgrnd, c.jpgrnd])) ...
    + ice_frac .* zso_capa;

% thermal conductivity of the layer calculated as the inverse of the sum of
% the inversed conductivities.
% thick_tt/k_tt_eff = thick_1/k_1 + thick_2/k_2 + thick_3/k_3
% Warning: It dos not include thermal exchange through air in pore and
% through water. In W/m/K.
zkappa = 1./ ...
    ( snow_frac_lay_1./zsn_condF(prhofirn, c) + ...
    snow_frac_lay_2./zsn_condF(prhofirn([2:c.jpgrnd, c.jpgrnd]), c) + ...
    ice_frac./zso_cond );

%% Calculate volumetric heat capacity and effective thermal conductivity
% by multiplying, resp. dividing, by each layer thickness to get to the 
% effective values

% calculating layer volume (without liquid water)
c.cdelV(1:c.jpgrnd) = psnowc(1:c.jpgrnd) *c .rho_water./prhofirn(1:c.jpgrnd) ...
    + psnic(1:c.jpgrnd)*c.rh2oice;
% !!! Why here not totalV ?

% Total heat capacity for each layer in W/K
zcapa_abs =zcapa .*c.cdelV ./ c.zdtime;

% The following used to go jk=1,c.jpgrnd-1
% Now we include c.jpgrnd, because it is needed below:

% Real distance between midpoints. Note repeated index for last layer.
cmidV =  (c.cdelV([2:c.jpgrnd, c.jpgrnd])+c.cdelV)./2;
% calculating effective thermal conductivity
zkappa_abs = zkappa ./ cmidV;
% !!! why her not *cmidV

% We have now calculated our new versions of zdz1, zdz2, zcapa and zkappa.
% Before introducing diffusion with sub-model layer, the old code was used 
% as it were. Now, we for some more:

% In the original version, this loop calculated c and d of c.jpgrnd-1 
% (which are in turn used above (in next time step) to prognose 
% tsoil(c.jpgrnd). Now we calculate c and d of c.jpgrnd.
% These are not used to prognose t of c.jpgrnd+1 (because this isn't 
% relevant) but to initialize the upwards calculation of c and d:

% Sublayer is all ice (zcapa = zso_capa) and of mass c.cdel(c.jpgrnd)
% giving physical thickness c.cdel(c.jpgrnd)*c.rh2oice:
% Update BV2017: using the temperature dependant heat capacity of ice
zcapa_abs_sublayer = zsn_capaF(c.rho_ice,pTdeep) * c.cdel(c.jpgrnd) ...
    * c.rho_water / c.rho_ice / c.zdtime;
% corresponds to zcapa_abs(c.jpgrnd+1)

z1 = zcapa_abs_sublayer + zkappa_abs(c.jpgrnd);

pgrndc(c.jpgrnd) = zcapa_abs_sublayer .* pTdeep ./ z1;
pgrndd(c.jpgrnd)=zkappa_abs(c.jpgrnd)./z1;
% PLA Tdeep (feb 2015) pTdeep changed in the above

% This loop went jk=c.jpgrnd-1,2,-1 (ie, calculating c and d for 3,2,1)
% Now, it goes   jk=c.jpgrnd,2,-1 (ie, calculating c and d for 4,3,2,1)
% It thus needs zdz1(c.jpgrnd) which is now also calculated above

for jk=c.jpgrnd:-1:2
    z1= 1/(zcapa_abs(jk) + zkappa_abs(jk-1) +  zkappa_abs(jk)*(1-pgrndd(jk)));

    pgrndc(jk-1)=(ptsoil(jk)*zcapa_abs(jk) +  zkappa_abs(jk)*pgrndc(jk))*z1;
    
    pgrndd(jk-1)=zkappa_abs(jk-1)*z1;
end
%   ---------------------------------------
%   COMPUTATION OFSIVE FLUX FROM GROUND AND
%   CALORIFIC CAPACITY OF THE GROUND:
%   ---------------------------------------------------------

pgrndhflx=zkappa_abs(1) * (pgrndc(1)          ...
    +(pgrndd(1)-1) * ptsoil(1));
pgrndcapc=(zcapa_abs(1) * c.zdtime+                ...
    c.zdtime * (1-pgrndd(1)) * zkappa_abs(1));
end

