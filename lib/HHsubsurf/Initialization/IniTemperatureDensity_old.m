function [snowthick, T_ice, rhofirn_ini,rho,snic_ini, snowc_ini, ...
    slwc_ini, z_icehorizon, Tsurf] = ...
    InitializationSubsurface(elev, T_obs, depth_thermistor, T_ice, ...
    Tsurf_obs,time, j, c)

% IniTemperatureDensity: Sets the initial state of the sub surface parameters:
% - snow depth
% - temperature profile
% - density profile
%
%
% Author: Dirk Van As (dva@geus.dk) & Robert S. Fausto (rsf@geus.dk)
% translated to matlab by Baptiste Vandecrux (bava@byg.dtu.dk)
% updated by Baptiste Vandecrux for integration of firn cores and
% temperature profiles as initial state of the subsurface column
%==========================================================================

snowthick=zeros(c.M,c.elev_bins);
%Initial snow depth
if c.elev_bins ==  1
    snowthick(1,j) = c.snowthick_ini ;
elseif elev(j) < c.ELA
    snowthick(1,j) = c.snowthick_ini + c.gradsnowthick*(elev(j)-c.ELA);
end

%% ========== Initial density profile ==================================
% Here we read the density profile from a core defined in a csv file in
% Input folder, average it on the weq depth scale used by the model, and
% give it as initial density profile for the model.  Below the depth of the
% given core, it is assumed to find ice.
%
% See script "InitialDensityProfileGenerator.m" in lib for creation of the
% csv file.
switch c.station
    case 'KAN-U'
%         if c.year(1) == 2012
            filename = './Input/Initial state/DensityProfile_KAN-U_2012.csv';
%         else
%             filename = './Input/Initial state/DensityProfile_KAN-U_1989.csv';
%         end
   case 'Miege'
        filename = './Input/Initial state/DensityProfile_Miege_2013.csv';
    case 'CP1'
        filename = './Input/Initial state/DensityProfile_CP1_1998.csv';
    case 'NUK_K'
        filename = './Input/Initial state/DensityProfile_NUK_K.csv';
    case 'Summit'
        filename = './Input/Initial state/DensityProfile_Summit_1990.csv';
    case 'NASA-SE'
        filename = './Input/Initial state/DensityProfile_NASA-SE_1998.csv';
    case 'NASA-E'
        filename = './Input/Initial state/DensityProfile_NASA-SE_1998.csv';  
    case 'NASA-U'
        filename = './Input/Initial state/DensityProfile_NASA-SE_1998.csv';
    case 'TUNU-N'
        filename = './Input/Initial state/DensityProfile_NASA-SE_1998.csv';  
    case {'DYE-2','DYE-2_long'}
        filename = './Input/Initial state/DensityProfile_DYE-2_1998.csv';
    case 'DYE-2_HQ'
        filename = './Input/Initial state/DensityProfile_DYE-2_HQ_2016.csv';
    otherwise
        disp('Missing initial density profile for requested station.');
end

out = dlmread(filename);
depth = out(:,1)/100; %real depth scale
olddensity = out(:,2); %density in kg/m^3 on the old real depth scale
if size(out,2)>2
    oldstrat = out(:,3); %density in kg/m^3 on the old real depth scale
else
    oldstrat=zeros(size(olddensity));
end

%the depthweq of each layer depend on the thickness_weq of all overlying
%layers. It is therefore necessary to fill missing data. Missing density
%are thus replaced by 350 kg/m^3.
olddensity(isnan(olddensity))=350;

thickweq = zeros(size(depth)); %thickness of each old layer in water eq
%calculating weq thickness of each layer
thickweq(1:length(depth)-1) = (depth(2:length(depth))-depth(1:length(depth)-1)) ...
    .*olddensity(1:length(depth)-1)/c.rho_water;
thickweq(end) = (depth(end)-depth(end-1))*olddensity(end)/c.rho_water;

oldscale_weq = zeros(size(depth)); %depthweq of each layer in original scale
% calculating the old weq depth scale
oldscale_weq(:) = cumsum(thickweq(:));

% calculates the new depth scale in mweq
depth_weq = zeros(size(c.cdel));
depth_weq = cumsum(c.cdel);
rhofirn_ini = zeros(size(c.cdel));
snowc_ini = zeros(size(c.cdel));
snic_ini = zeros(size(c.cdel));
slwc_ini = zeros(size(c.cdel));

% we limit the new scale to the values that are given within the oldscale
newscale_weq = depth_weq(depth_weq<oldscale_weq(end));
olddensity = olddensity(oldscale_weq<=newscale_weq(end));
oldstrat = oldstrat(oldscale_weq<=newscale_weq(end));
oldscale_weq = oldscale_weq(oldscale_weq<=newscale_weq(end));
newscale_weq = depth_weq(depth_weq<oldscale_weq(end));

% Since the oldscale is denser than the new one, we need to make the union
% of both scales and then average the density for each section of the new
% scale.
mergedscale = union(oldscale_weq,newscale_weq);
density_mergedscale = interp1(oldscale_weq,olddensity,mergedscale,'nearest');
ice_perc_mergedscale = interp1(oldscale_weq,oldstrat,mergedscale,'nearest');

% In the coming section we further process the core to increase resolution.
% The stratigraphy is used to know the percentage of ice within each core
% section. Then it is used to initiate snic. The firn remaining in each
% section is used to initiate snowc and the weight of the section (once
% subtracted the weight of ice content) is used to initiate prhofirn.

    % calculating the thickness of each layer in merged scale
    thick_mergedscale_weq = mergedscale;
    thick_mergedscale_weq(2:length(mergedscale)) = mergedscale(2:end) - mergedscale(1:end-1);
    thick_mergedscale_m = thick_mergedscale_weq * c.rho_water ...
        ./ density_mergedscale;

    % calculating the volume and mass of ice contained in each section
    icevol_mergedscale_m = thick_mergedscale_m .* ice_perc_mergedscale/100;
    icemass_mergedscale_kg = icevol_mergedscale_m * 800; 
    % in line above 600 is taken for density of ice, it is low but it was
    % noticed that putting higher values  leads to very low density for the
    % firn left in the section.

    % calculating the volume and mass of firn left in the section excluding
    % the ice content
    firnvol_mergedscale_m = thick_mergedscale_m .* (1-ice_perc_mergedscale/100);
    firnmass_mergedscale_kg = density_mergedscale .* thick_mergedscale_m - icemass_mergedscale_kg;
    firndensity_mergedscale = firnmass_mergedscale_kg ./ firnvol_mergedscale_m;
    firndensity_mergedscale(firnvol_mergedscale_m==0)=500; % giving an aritificial density of 500 to the firn in sections only made of ice
    
    ind_bins = discretize(mergedscale,[0; depth_weq]);
    nbins = length(unique(ind_bins));
    mean_firndensity = zeros(nbins, 1);

    % now we either sum or average the ice content, snow content and firn
    % density on the model depth scale
    for ii = 1:nbins
        ind_in_bin = (ind_bins == ii);
%         if sum(firnvol_mergedscale_m(ind_in_bin)) ~= 0
%             mean_firndensity(ii) = sum(firndensity_mergedscale(ind_in_bin) ...
%                 .* firnvol_mergedscale_m(ind_in_bin)) ...
%                 / sum(firnvol_mergedscale_m(ind_in_bin));
%         else
%              mean_firndensity(ii) = 500;
%         end
    mean_firndensity(ii) = sum(density_mergedscale(ind_in_bin).*thick_mergedscale_m(ind_in_bin)) ...
        /sum(thick_mergedscale_m(ind_in_bin));

%         snic_ini(ii) = min(sum(icemass_mergedscale_kg(ind_in_bin))./c.rho_water,c.cdel(ii));
%         snowc_ini(ii) = c.cdel(ii)-snic_ini(ii);
        snowc_ini = c.cdel;
        snic_ini(ii) = 0;
        slwc_ini(ii) = 0;
        if snowc_ini(ii)<0
            eweojjnjv;
        end
    end
 
% giving the observed density (on appropriate scale) as initial value
% for the subsurface density profile
rhofirn_ini(1:length(mean_firndensity)) = min(max(300, mean_firndensity),900);

% fills the rest of the density profile with ice
if strcmp(c.station,'NUK_K')
        x =[depth_weq(1:length(mean_firndensity)); 50; 60];
        y = [mean_firndensity; 900;900];
        fo = fit(x, y, 'poly2');

        rhofirn_ini(length(mean_firndensity)+1:length(c.cdel),1) = ...
            860;
        rhofirn_ini(rhofirn_ini>c.rho_ice) = c.rho_ice;
        snic_ini(length(mean_firndensity)+1:length(c.cdel),1) = c.cdel(length(mean_firndensity)+1:length(c.cdel),1);
        snowc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0;
        slwc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0;
else
    if length(mean_firndensity)<length(c.cdel)
        x =[depth_weq(1:length(mean_firndensity)); 50; 60];
        y = [mean_firndensity; 900;900];
        fo = fit(x, y, 'poly2');

        rhofirn_ini(length(mean_firndensity)+1:length(c.cdel),1) = ...
            fo(depth_weq(length(mean_firndensity)+1:length(c.cdel)));
        rhofirn_ini(rhofirn_ini>c.rho_ice) = c.rho_ice;
        snic_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0;
        snowc_ini(length(mean_firndensity)+1:length(c.cdel),1) = c.cdel(length(mean_firndensity)+1:length(c.cdel),1) ;
        slwc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0;
    end
end
% figure
% plot(depth_weq, rhofirn_ini)

rho =  (c.cdel *c.rho_water)./...
    ((snowc_ini *c.rho_water)./rhofirn_ini + (snic_ini *c.rho_water)./c.rho_ice);

%Snow to ice transition
z_icehorizon = min(c.z_ice_max, floor(snowthick(1,j)/c.dz_ice));

%Overwrite density in the snow pack with a given value for snow density
% -> here not needed because snow density should be given in initial
% density profile

% if snowthick(1,j) > 0
%     rho(1:z_icehorizon,1) = c.rho_snow(1,j);
%     %+ (c.rho_ice-c.rho_snow(j,0))*z*c.dz_ice/snowthick(j,0)   % including snow in the density array
% end

% Removing densities greater than ice densities (just in case)
toodense = (rhofirn_ini > c.rho_ice);
if sum(toodense) >0
    rhofirn_ini(toodense) = c.rho_ice;
end

if c.verbose == 1
    f = figure('Visible','on');
    subplot(1,2,1)
    plot(oldscale_weq,olddensity,'--')
    hold on
    stairs([0; depth_weq(1:end-1)],rhofirn_ini, 'LineWidth',2)
    stairs([0; depth_weq(1:end-1)],rho, 'LineWidth',2)
    plot([0, depth_weq(end-1)],[c.rho_ice, c.rho_ice],'k')
    legend('data','firn density','bulk density','Location','South')
    xlabel('Depth (m weq)')
    ylabel('Density (kg/m^3)')
    title('Initial density profile')
    view([90 90])
end
%% ========== initial surface temperature ================
% Initial Tsurf: reads initial surface temperature from AWS
Tsurf_ini = Tsurf_obs; %in deg K
%could also be defined as user-defined param in c struct

Tsurf = (Tsurf_ini) * ones(c.M,c.elev_bins);

%% ================= Initial temperature profile =========================
% same as density, we give initial temperature to the subsurface according
% to a csv file located in the Input folder
%
% See script "InitialTemperatureProfileGenerator.m" in lib folder for creation of
% the csv file.
   
% if there is thermistor record for the first time step, then reads 
% initial subsurface conditions from AWS data
oldtemp = [];
time_dt = datenum(time,1,1);

for i = 1:24*7 %we look at the week following the installation
    if sum(~isnan(T_obs(i,:)))>1
        depth= depth_thermistor(i,~isnan(T_obs(1,:)))'; %in m
        depth = depth(depth~=0);
        [depth,ind_sorted] = sort(depth);
        
        oldtemp = T_obs(i,~isnan(T_obs(1,:)))';
        oldtemp = oldtemp(depth~=0);
        oldtemp = oldtemp(ind_sorted);
        date_Tstring = datestr(time_dt(i));
        break
    end
end
if isempty(oldtemp) || strcmp(c.station,'NUK_K')    
    %reads initial subsurface conditions from file
    out = dlmread('./Input/Initial state/InitialTemperatureProfile.csv');
    switch c.station
%         case 'Miege'
%         out = dlmread('./Input/Initial state/InitialTemperatureProfile_Miege.csv');
        case 'NUK_K'
        out = dlmread('./Input/Initial state/InitialTemperatureProfile_NUK_K.csv');

    end

    date_Tstring = 'manualy chosen profile';
    depth = out(:,1); %old depth scale in m
    oldtemp = out(:,2); %density in kg/m^3 on the old real depth scale
    
end
oldtemp(depth<0 ) =[];
depth(depth<0 ) =[];

% Preparing the model's depth scale
    % calculates the new depth scale in mweq
    depth_weq = zeros(size(c.cdel));
    depth_weq(:) = cumsum(c.cdel(:));

    % calculates the new depth scale in real m
    depth_act = cumsum(c.cdel .*c.rho_water ./rhofirn_ini(:,1));

    % Here we add an initial surface temperature
    depth_act = [0; depth_act];
    depth_weq = [0; depth_weq];
    
% Prepraring observation's depth scale
    if depth(1) ~= 0
        %if the input density profile does not contain surface temperature,
        %then we use Tsurf_in to initiate the temperature profile
        depth = [0; depth];
        oldtemp = [Tsurf_ini - c.T_0; oldtemp];
    end
    oldtemp(depth > depth_act(end)) = [];
    depth(depth > depth_act(end))=[];

    % the old scale is converted from m to mweq by interpolation
    oldscale_weq = interp1(depth_act,depth_weq,depth);

    % we limit the new scale to the values that are given within the oldscale
    newscale_weq = depth_weq(depth_weq<=oldscale_weq(end));

    % the temperature is interpolated on each stamp of the new scale
    newtemp = interp1(oldscale_weq,oldtemp,newscale_weq);
    newtemp = newtemp + c.T_0; %going back to K

    % giving the observed temperature (on appropriate scale) as initial value
    % for the subsurface temperature profile
    % There might be a shift to apply depending on whether the first value in
    % subsurface column represent the temp at depth 0 (=Tsurf) or at depth 
    % c.cdel(1). To be asked.
    T_ice(1:length(newtemp),1,j) = newtemp;

    if length(newtemp)<length(c.cdel)
        d1 = length(newtemp);
        T1 = newtemp(d1);
        d2 = c.jpgrnd;
        T2 = c.Tdeep_AWS + c.T_0;

        ind = length(newtemp)+1:(c.jpgrnd-1);
        T_ice(ind,1,j) = ...
            (T1-T2)/(d2-d1)^2*(d2-ind).^2 + T2;
        T_ice(c.jpgrnd,1,j) = c.Tdeep_AWS+ c.T_0;
    end
if c.verbose == 1
    subplot(1,2,2)
    scatter(oldscale_weq,oldtemp,'o')
    hold on
    stairs(depth_weq(1:end-1),T_ice(:,1,j)-c.T_0, 'LineWidth',2)
    legend(sprintf('data from %s',date_Tstring),'after interpolation','Location','South')
    xlabel('Depth (m weq)')
    ylabel('Temperature (deg C)')
    title(sprintf('Initial temperature profile\n on %s',datestr(time_dt(1))))
    view([90 90])
    print(f, sprintf('%s/IniTempDens',c.OutputFolder), '-dpng');
disp('Initial subsurface profiles generated')

end
% removing non-freezing temperatures (just in case)
subsurfmelt = find(T_ice(:,1,j) > c.T_0);
if sum(subsurfmelt )>0
    T_ice(subsurfmelt,1,j) = c.T_0;
end

end
