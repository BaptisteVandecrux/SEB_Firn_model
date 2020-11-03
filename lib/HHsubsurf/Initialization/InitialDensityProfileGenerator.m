% This script helps to generate a density profile to be used as initial
% state for the subsurface model
%
% Author: Baptiste Vandecrux (bava@byg.dtu.dk)
% ========================================================================

% clear all
% close all
load ../matlab_functions/Core_all


station = 'EGP';
% Plot in the site you want
%  PlotCore(Core,'Site','Camp Century')

% Find the index of the core you want to use
% CoreList(CoreAvg)
% search by name rather than by index since index can change in the future

switch station
    case 'CP1'
        ind = FindCore(Core,'Name','CORE 6945');
    case 'DYE-2'
        ind = FindCore(Core,'Name','DYE2 1998 core B');
    case 'DYE-2_HQ'
        ind = FindCore(Core,'Name','core_10_2016');
    case 'KAN-U'
        ind = FindCore(Core,'Name','core_1_2012');
    case 'NASA-SE'
        ind = FindCore(Core,'Name','CORE 6642 (B)');
    case 'EGP'
        ind = FindCore(Core,'Name','NEGIS');
%        ind_new = length(Core)+1;
%        Core{ind_new} = Core{ind1};
%        Core{ind_new}.Info.Name = 'NASA-SE_bapt';
%        Core{ind_new}.Data.Density(50:110) = NaN;
%        ind_nan = isnan(Core{ind_new}.Data.Density);
%        Core{ind_new}.Data.Density(ind_nan) = ...
%            interp1(Core{ind_new}.Data.Depth(~ind_nan),...
%            Core{ind_new}.Data.Density(~ind_nan),...
%            Core{ind_new}.Data.Depth(ind_nan));
%        figure
%        OverlapPlot(Core,[ind1 ind_new])
%         ind = ind_new;
        
%         ind2=FindCore(Core,'NearestCodeLocation','NASA-SE');
%         figure
%         PlotCore(Core,'CoreNumber',ind2)
%        figure
%        OverlapPlot(Core,[114 23])
    case 'Summit'
        ind1 = FindCore(Core,'Name','T99_1990');
       ind2 = FindCore(Core,'Name','Grip1991Shallow');
       ind_new = length(Core)+1;
       Core{ind_new} = Core{ind2};
       Core{ind_new}.Info.Name = 'T99+GRIP_1990';
       Core{ind_new}.Data.Density(1:length(Core{ind1}.Data.Density)) = ...
           Core{ind1}.Data.Density;
       ind_nan = isnan(Core{ind_new}.Data.Density);
       Core{ind_new}.Data.Density(ind_nan) = ...
           interp1(Core{ind_new}.Data.Depth(~ind_nan),...
           Core{ind_new}.Data.Density(~ind_nan),Core{ind_new}.Data.Depth(ind_nan));
       
        ind = ind_new;
    case 'Miege'
        ind = 120; %or 8
    case 'NASA-U'
        ind = FindCore(Core,'Name','CORE 7347');
%         ind2 = FindCore(Core,'Name','NASA-U_Henderson');

    case 'SouthDome'
        ind = FindCore(Core,'Name','S. Dome Core A');
%         ind = FindCore(Core,'Name','S. Dome Core B');
    case 'NASA-E'
        ind1 = FindCore(Core,'Name','NASA East Core A');
        ind2 = FindCore(Core,'Name','NASA East Core B');
        ind_new = length(Core) +1;
        Core{ind_new} = Core{ind1};
        Core{ind_new}.Info.Name = 'NASA East combined';
        tmp = length(Core{ind2}.Data.Density);
        Core{ind_new}.Data.Density(1:tmp) = Core{ind2}.Data.Density;

    case 'Saddle'
        ind = FindCore(Core,'Name','N. Dye 3 (Saddle) - B');         
        ind2 = FindCore(Core,'Name','N. Dye 3 (Saddle) - A');
        figure
        OverlapPlot(Core,[ind ind2]);
    case 'TUNU-N'
        %%
%         ind2 = FindCore(Core,'Name','B19_NGT19_1994');
%         ind2 = FindCore(Core,'Name','Tunu-S7.5');
%         tunus7.5 lower than tunu1
%         B19 lower than tunu1
%         tunuN25 lower than tunu1
%         tunuE25 lower than tunu1
%         tunuW25 lower than tunu1


        ind1 = FindCore(Core,'Name','Tunu-W25');
        ind2 = FindCore(Core,'Name','Tunu-E50');

        ind = FindCore(Core,'Name','Tunu-1');

        figure
        OverlapPlot(Core,[ind ind2]);
        %%
    case 'NGRIP'
        ind = FindCore(Core,'Name','NG97S2~1-3bag');
%         ind = FindCore(Core,'Name','NGRIP2001S5');
    case 'GITS'
        ind = FindCore(Core,'Name','Camp Century');
%         ind = FindCore(Core,'Name','NGRIP2001S5');
       
end

% Creating density profile and writing it into the Input folder
depth = Core{ind}.Data.Depth/100;
density = Core{ind}.Data.Density;
ice_perc = Core{ind}.Data.Type_perc;

% if strcmp(station,'KAN-U')
%     depth = depth';
%     density = density';
%     ice_perc = ice_perc';
%     ice_perc(length(ice_perc):length(density)) = 0;
% end
% DensProfile = [depth, density, ice_perc];
filename = sprintf('./Input/Initial state/DensityProfile_%s_%i.csv',station,Core{ind}.Info.DateCored.Year);
    M  = [depth', density'];
    M_table = array2table(M,'VariableName', {'depth_m', 'density_kgm3'});

    writetable(M_table,filename,'Delimiter',';')

fprintf('Initial density profile was generated from core %s and placed in Input folder.\n',Core{8}.Info.Name)
PlotCore(Core,'CoreNumber',ind);





