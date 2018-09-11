function  [H_comp, GF, GFsurf, GFsubsurf, H_melt_weq,H_rain,...
H_subl, H_surf,H_surf_weq, H_snow, L, LHF, meltflux, meltflux_internal, ...
rho, rho_snow, runoff,...
SHF, SRnet,SRout_mdl,  T_ice, grndc, grndd , slwc, snic, snowc, zrogl, ...
pdgrain, refreezing, snowbkt, theta_2m,q_2m,ws_10m, Tsurf, snowthick,...
z_icehorizon,Re,Ri,err]...
= IniVar(c)
%IniVar: Initiate all modelled variable.
%Matlab does run faster when memory is preallocated.
%
% Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================

% DECLARATION OF ARRAYS OF MODELLED VARIABLES -----------------------------------------------------------------------
GF         = zeros(c.jpgrnd);
GFsurf        = zeros(c.M,c.elev_bins);
GFsubsurf        = zeros(c.M,c.elev_bins);
H_melt_weq        = zeros(c.M,c.elev_bins);
H_rain        = zeros(c.M,c.elev_bins);
H_subl        = zeros(c.M,c.elev_bins);
H_surf        = zeros(c.M,c.elev_bins);
H_surf_weq        = zeros(c.M,c.elev_bins);
H_snow        = zeros(c.M,c.elev_bins);
H_comp       = zeros(c.M,c.elev_bins);
L         = zeros(c.M,c.elev_bins);
LHF         = zeros(c.M,c.elev_bins);
meltflux      = zeros(c.M,c.elev_bins);
meltflux_internal = zeros(c.M,c.elev_bins);
refreezing      = zeros(c.jpgrnd,c.M,c.elev_bins);
rho         = zeros(c.jpgrnd,c.M);
rho_snow      = zeros(c.M,c.elev_bins);
runoff        = zeros(c.M,c.elev_bins);
SHF         = zeros(c.M,c.elev_bins);
theta_2m         = zeros(c.M,c.elev_bins);
Tsurf         = zeros(c.M,c.elev_bins);
q_2m         = zeros(c.M,c.elev_bins);
Re         = zeros(c.M,c.elev_bins);
Ri         = zeros(c.M,c.elev_bins);
err         = zeros(c.M,c.elev_bins);
ws_10m         = zeros(c.M,c.elev_bins);
SRnet       = zeros(c.jpgrnd,1);
SRout_mdl = zeros(c.jpgrnd,1);
T_ice       = zeros(c.jpgrnd,c.M,c.elev_bins);   % sub-surface temperatures (snow and ice)
    snowthick=zeros(c.M,c.elev_bins);
z_icehorizon = NaN;
grndc      =zeros(c.jpgrnd,1);
grndd      =zeros(c.jpgrnd,1);
slwc       =zeros(c.jpgrnd,1);
snic       =zeros(c.jpgrnd,1);
snowc      = c.cdel;
pdgrain      =zeros(c.jpgrnd,1);
% DECLARATION OF ARRAYS OF MODELLED VARIABLES -----------------------------------------------------------------------
zrogl=0;
snowbkt = 0;
end
