% Main script for running the surface - subsurface model
% Here you can define which year you want to compute, define your
% parameters, plot output variables, evaluate model performance....
%
% Author: Baptiste Vandecrux (bava@byg.dtu.dk)
% ========================================================================
clearvars
close all
clc
addpath(genpath('lib'))
addpath(genpath('Input'),genpath('Output'))

% All parameters are defined in a csv files in Input folder, however it is
% possible to overwrite the default value by defining them again in the
% "param" struct hereunder.

% Increasing resolution towards the sufrace
% param.cdel = zeros(32,1);
% param.cdel = 3.^([1:32]'-29)+0.15;

% High resolution grid (comment if not needed)
NumLayer = 200;
param.z_max = 40;
param.dz_ice = param.z_max/NumLayer;
param.verbose = 1;
param.lim_new_lay = 0.04;

% Heterogeneous precolation (comment if not needed)
% param.hetero_percol = 0;
% param.hetero_percol_p = 0.01;
% param.hetero_percol_frac = 0.5;

param.snowthick_ini= 0.1* 350 / 1000; %Initial snow thickness in m weq

param.calc_CLliq = 1; % uses Coleou and Lesaffre to calculate irreduscible water content

param.do_no_darcy = 0; % = 0 for using the new darcy-like bucket scheme
param.whwice = 0.1;  % ratio between ice layer horizontal span and horizontal spacing from other ice layers
param.Ck  =  1; % refreezing speed : fraction of the potential 
               % refrozen mass actually refrozen per time step
% param.liqmax          =      0.2;

param.rho_snow_scheme = 0;
% choice for fresh snow density see function IniRhoSnow
% 0 -> constant = 315; % Greenland-wide mean according to Fausto et al 201
% 1 -> dependant on mean temperature (Reeh et al. 2005)
% 2 -> dependant on mean temperature (Kuipers-Munneke et al. 2015)
% 3 -> dependant on elevation latitude and longitude (RSF unpub)
% 4 -> dependant on hourly temperature and wind speed (Liston et al., 2007)

param.precip_scheme = 3;
% choice for precipitation scheme see function Precipitation
% 1 -> when LRin_AWS >= c.sigma*T_AWS.^4, constant rate c.precip_rate (Dirk, unpub.)
param.prec_rate       =      0.001;
% 2 -> when RH > 80%, with intensity scaled by (RH-80)/20, baseline rate
% c.precip_rate for 100% humidity (Liston and Sturm 1998)
% 3 -> Using surface height measurements (smoothed over 1 week) to quantify
% the snowfall

param.ConductionModel = 0;      % if 1 does CONDUCTION ONLY
% In the conduction model, the temperature profile is reseted every night
% at 2am local time usingg thermistor string readings

param.year    = 0;
% by defining param.year, the model will be run only for that melt year
% (i.e. 1st. april to 1st april) however you still need to make sure that 
% "rows"  is set so that the appropriate values will be read in the weather 
% data. To run the model only from 1 to rows, just leave param.year=0.

% param.perturb_rho_ini = -40; % this option allows to perturbate the initial 
                            % density profile this shift will be applied on
                            % the upper 20 m of cores.

%    param.a_dens = 30.25;
%    param.b_dens = 0.7;
param.track_density = 1;
param.avoid_runoff = 1;
param.THF_calc = 1;

station_list = {'DYE-2_long','CP1', 'Summit','NASA-SE'};
% 
% 
for i =1:length(station_list)
param.station =  station_list{i}; %'NASA-SE';

switch param.station
    case 'CP1'
        param.InputAWSFile = 'data_CP1_1998-2010.txt';
    case 'DYE-2'
        param.InputAWSFile = 'data_DYE-2_1998-2015.txt';
%         param.InputAWSFile = 'data_DYE-2_combined_hour.txt';
%         param.InputAWSFile = 'data_DYE-2_combined_hour_cor.txt';
%         param.InputAWSFile = 'data_DYE-2_restricted_hour.txt';
    case 'DYE-2_long'
        param.InputAWSFile = 'data_DYE-2_1998-2015.txt';

    case 'Summit'
%         param.InputAWSFile = 'data_Summit_1990-2015.txt';
        param.InputAWSFile = 'data_Summit_2000-2015.txt';
%         param.InputAWSFile = 'data_Summit_combined_hour.txt';
        
    case 'NASA-SE'
        param.InputAWSFile = 'data_NASA-SE_1998-2015.txt';
      otherwise
        disp('Missing data file for the requested station.');
end

% When you add sites
% 1) in Main.m: define path of InputAWSFile
% 2) in IniTemperatureDensity.m: Define path for density profile
% 3) Make sure all information is reported in Input/Constants/InfoStation.csv file

%%  ========= Run model with name tag of your choice =======================
[RunName, c] = HHsubsurf(param);
end
% 