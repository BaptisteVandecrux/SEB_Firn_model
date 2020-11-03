function [T_ice, rhofirn_ini,rho,snic_ini, snowc_ini, ...
    slwc_ini, graind_out] = ...
    InitializationSubsurface(T_obs, depth_thermistor, T_ice, ...
    time, Tsurf_ini, j, c)

% InitializationSubsurface: Sets the initial state of the sub surface parameters:
% - snow depth
% - temperature profile
% - density profile
%
%
% Author: Baptiste Vandecrux (bava@byg.dtu.dk)
%==========================================================================


%% ========== Initial density profile ==================================
% Here we read the density profile from a core defined in a csv file in
% Input folder, average it on the weq depth scale used by the model, and
% give it as initial density profile for the model.  Below the depth of the
% given core, it is assumed to find ice.
%
% See script "InitialDensityProfileGenerator.m" in lib for creation of the
% csv file.
if c.retmip == 0
    switch c.station
        case {'DYE-2','DYE-2_long'}
            filename = './Input/Initial state/density/RetMIP_density_Dye-2 1998.csv';
        case {'CP1', 'Summit', 'NASA-SE', 'NASA-E', 'NASA-U', 'SouthDome', 'Saddle', 'TUNU-N'}
            filename = ['./Input/Initial state/density/DensityProfile_', c.station,'_1998.csv'];
%             filename = './Input/Initial state/density/RetMIP_density_Summit 1990.csv';
        case 'NGRIP'
            filename = './Input/Initial state/density/DensityProfile_NGRIP_1997.csv';  
        case 'GITS'
            filename = './Input/Initial state/density/DensityProfile_GITS_1963.csv';  
        case 'DYE-2_HQ'
            filename = './Input/Initial state/density/RetMIP_density_Dye-2 2016.csv';
        case 'NUK_K'
            filename = './Input/Initial state/density/DensityProfile_NUK_K.csv';
        case 'EGP'
            filename = './Input/Initial state/density/DensityProfile_EGP_2012.csv';
        case 'KAN-U'
    %         if c.year(1) == 2012
                filename = './Input/Initial state/density/RetMIP_density_KAN-U 2012.csv';
    %         else
    %             filename = './Input/Initial state/density/DensityProfile_KAN-U_1989.csv';
    %         end
       case 'Miege'
            filename = './Input/Initial state/density/RetMIP_density_Firn Aquifer.csv';
        otherwise
            if c.elev_AWS < 1400
                filename = './Input/Initial state/density/DensityProfile_NUK_K.csv';
            elseif c.elev_AWS < 1800
                filename = './Input/Initial state/density/RetMIP_density_Dye-2 2016.csv';
            else
                filename = './Input/Initial state/density/DensityProfile_Summit_1998.csv';
            end
            warning('Missing initial density profile for requested station.');
    end
else
    filename = ['.\RetMIP\Input files\density\RetMIP_density_', c.station,'.csv'];
end
disp('File used for firn density initialization:')
disp(filename)

delimiter = ';';
startRow = 1;
formatSpec = '%f%f%[^\n\r]';

fileID = fopen(filename,'r');
try dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
    'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
catch me
    startRow = 2;
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
        'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
end
fclose(fileID);

out = [dataArray{1:end-1}];

depth = out(:,1); %real depth scale
olddensity = out(:,2); %density in kg/m^3 on the old real depth scale
if size(out,2)>2
    oldstrat = out(:,3); %density in kg/m^3 on the old real depth scale
else
    oldstrat=zeros(size(olddensity));
end

if c.perturb_rho_ini ~= 0
    old_thickness = depth;
    old_thickness(2:end ) = depth(2:end) - depth(1:end-1);
    [~, ind_20] = min(abs(depth-20));

    olddensity_avg = sum(olddensity(1:ind_20).*old_thickness(1:ind_20)) / depth(ind_20);
    perturb_olddensity = olddensity(1:ind_20)...
        + c.perturb_rho_ini * depth(ind_20)/length(1:ind_20)./ old_thickness(1:ind_20);
    perturb_olddensity_avg = sum(perturb_olddensity.*old_thickness(1:ind_20))...
        / depth(ind_20);
    perturb_olddensity(perturb_olddensity>917)=917;

    f = figure('Visible','on','OuterPosition',[0 0 8 18]);
    plot(olddensity,-depth)
    hold on
    plot(perturb_olddensity,-depth(1:ind_20))
    ylabel('Depth (m)')
        xlabel('Density (kg m^{-3})','Interpreter','tex')
    legend('Original','After perturbation')
    title(sprintf('Initial density perturbated\non average by %0.1f kg/m-3 \nin the upper 20 m.',perturb_olddensity_avg - olddensity_avg))

    olddensity(1:ind_20) = perturb_olddensity;
end
%the depthweq of each layer depend on the thickness_weq of all overlying
%layers. It is therefore necessary to fill missing data. Missing density
%are thus replaced by 350 kg/m^3.
olddensity(isnan(olddensity))=350;

thickweq = zeros(size(depth)); %thickness of each old layer in water eq
%calculating weq thickness of each layer
% tick_old_act =  depth;
% tick_old_act(2:end)  = depth(2:end) - depth(1:end-1);
% thickweq = thick_old_act .* olddensity /c.rho_water;
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
depth = depth(oldscale_weq<=newscale_weq(end));
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
    
    ind_bins = discretize(mergedscale,[0; depth_weq],'IncludedEdge','right');
    nbins = length(unique(ind_bins));
    mean_firndensity = zeros(nbins, 1);

    % now we either sum or average the ice content, snow content and firn
    % density on the model depth scale
    for ii = 1:nbins
        ind_in_bin = (ind_bins == ii);

    mean_firndensity(ii) = sum(density_mergedscale(ind_in_bin).*thick_mergedscale_m(ind_in_bin)) ...
        /sum(thick_mergedscale_m(ind_in_bin));

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
        rhofirn_ini(length(mean_firndensity)+1:length(c.cdel),1) = ...
            860;
        rhofirn_ini(rhofirn_ini>c.rho_ice) = c.rho_ice;
        snic_ini(length(mean_firndensity)+1:length(c.cdel),1) = c.cdel(length(mean_firndensity)+1:length(c.cdel),1);
        snowc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0;
        slwc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0;
else
    if length(mean_firndensity)<length(c.cdel)
        x =[depth_weq(1:length(mean_firndensity)); 30; 70];
        y = [mean_firndensity; 830;830];
        p = polyfit(x,y,2);
        fo = @(x) p(1)*x.^2 + p(2)*x + p(3);
        
        rhofirn_ini(length(mean_firndensity)+1:length(c.cdel),1) = ...
            fo(depth_weq(length(mean_firndensity)+1:length(c.cdel)));
        rhofirn_ini(rhofirn_ini>c.rho_ice) = c.rho_ice;
        snic_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0;
        snowc_ini(length(mean_firndensity)+1:length(c.cdel),1) = c.cdel(length(mean_firndensity)+1:length(c.cdel),1) ;
        slwc_ini(length(mean_firndensity)+1:length(c.cdel),1) = 0;
    end
end

rho =  (c.cdel *c.rho_water)./...
    ((snowc_ini *c.rho_water)./rhofirn_ini + (snic_ini *c.rho_water)./c.rho_ice);

% Removing densities greater than ice densities (just in case)
toodense = (rhofirn_ini > c.rho_ice);
rhofirn_ini(toodense) = c.rho_ice;

% calculates the new depth scale in real m
depth_act = cumsum(c.cdel .*c.rho_water ./rhofirn_ini(:,1));

if c.verbose == 1
    f = figure('Visible','on','OuterPosition',[0 0 8 18]);
stairs([0; depth(1:end-1)],olddensity,'LineWidth',2)
    hold on
    stairs([0; depth_act(1:end-1)],rhofirn_ini, 'LineWidth',1.5)
    plot([0, depth_weq(end-1)],[c.rho_ice, c.rho_ice],'k')
    axis tight
    box on
    legend('Observations','Model initial state',...'bulk density',
        'Location','NorthOutside')
    legend boxoff
    xlabel('Depth (m)')
    ylabel('Density (kg m^{-3})','Interpreter','tex')
    view([90 90])
    title(c.station)
        print(f,sprintf('%s/Initial_rho.tif',c.OutputFolder),'-dtiff')

end

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
clearvars out

    %reads initial subsurface conditions from file
if c.retmip == 0
    switch c.station
        case 'DYE-2_HQ'
        filename = './Input/Initial state/temperature/RetMIP_temperature_Dye-2 HQ.csv';
        case {'DYE-2','DYE-2_long'}
        filename = './Input/Initial state/temperature/RetMIP_temperature_Dye-2.csv';
        case {'Summit', 'NASA-E', 'NGRIP'}
        filename = './Input/Initial state/temperature/RetMIP_temperature_Summit.csv';
        case {'NASA-U', 'KAN-U'}
        filename = ['./Input/Initial state/temperature/RetMIP_temperature_', c.station,'.csv'];
        case 'GITS'
            filename = './Input/Initial state/temperature/GITS_temperature.csv';  
        case 'NUK_K'
            filename = './Input/Initial state/temperature/InitialTemperatureProfile_NUK_K.csv';
        case 'Miege'
            filename = './Input/Initial state/temperature/RetMIP_temperature_Firn Aquifer.csv';

    
        otherwise
            disp('Using thermistor data')
            for i = 1:24*14 %we look at the week following the installation
                if sum(~isnan(T_obs(i,:)))>1
                    if sum(~isnan(T_obs(i,:)))<3
                        continue
                    end
                    if sum(~isnan(depth_thermistor(i,~isnan(T_obs(i,:)))'))<3
                        continue
                    end
                    depth= depth_thermistor(i,~isnan(T_obs(i,:)))'; %in m
                    depth = depth(depth~=0);
                    [depth,ind_sorted] = sort(depth);
                    out(:,1) = depth;
                    oldtemp = T_obs(i,~isnan(T_obs(i,:)))';
                    oldtemp = oldtemp(depth~=0);
                    oldtemp = oldtemp(ind_sorted);
                    out(:,2) = oldtemp;
                    date_Tstring = datestr(time_dt(i));
                    break
                end
            end
            filename = [];
        
        if isempty(oldtemp)
            disp('no thermistor data availble')
            filename = './Input/Initial state/temperature/InitialTemperatureProfile.csv';
        end
        
    end

else
    filename = ['.\RetMIP\Input files\temperature\RetMIP_temperature_', c.station,'.csv'];
end

    if ~isempty(filename)
        disp('File used for firn temperature initialization:')
        disp(filename)
                delimiter = ';';
        startRow = 2;
        formatSpec = '%f%f%[^\n\r]';

        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
            'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        fclose(fileID);

        out = [dataArray{1:end-1}];
    end

% date_Tstring = 'manualy chosen profile';
depth = out(:,1); %old depth scale in m
oldtemp = out(:,2); %temperature on the old real depth scale


oldtemp(depth<0 ) =[];
depth(depth<0 ) =[];
  depth(isnan(oldtemp))=[];
  oldtemp(isnan(oldtemp))=[];
  
% Prepraring observation's depth scale
    if depth(1) ~= 0
        %if the input density profile does not contain surface temperature,
        %then we use Tsurf_in to initiate the temperature profile
        depth = [0; depth];
        oldtemp = [Tsurf_ini - c.T_0; oldtemp];
    end
    if depth_act(end)>depth(end)
        depth = [depth; depth_act(end)];
        oldtemp = [oldtemp; oldtemp(end)];
    end
        
    [newtemp] = ConvertToGivenDepthScale(depth, oldtemp, depth_act,'linear');

    newtemp(isnan(newtemp))=[];
    newtemp = newtemp + c.T_0; %going back to K

    % giving the observed temperature (on appropriate scale) as initial value
    % for the subsurface temperature profile
    % There might be a shift to apply depending on whether the first value in
    % subsurface column represent the temp at depth 0 (=Tsurf) or at depth 
    % c.cdel(1). To be asked.
    T_ice(1:length(newtemp),1,j) = newtemp;

    % EXTRAPOLATING UNTIL DEEP FIRN TEMPERATURE
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
    % removing non-freezing temperatures (just in case)
    subsurfmelt = find(T_ice(:,1,j) > c.T_0);
    if sum(subsurfmelt )>0
        T_ice(subsurfmelt,1,j) = c.T_0;
    end
    
if c.verbose == 1
    f = figure('Visible','on','OuterPosition',[0 0 8 18]);
    scatter(depth,oldtemp,'o','fill')
    hold on
    stairs(depth_act, T_ice(:,1,j)-c.T_0, 'LineWidth',2)
    scatter(depth(end),oldtemp(end),'s','LineWidth',6)
    legend('Observations','Model initial state',...
        'Prescribed deep \newlinefirn temperature',...
        'Interpreter','tex',    'Location','NorthOutside')

    ylabel('Temperature (^oC)','Interpreter','tex')
    legend boxoff
    axis tight
    box on
    view([90 90])  
    title(c.station)  
    print(f,sprintf('%s/Initial_temp',c.OutputFolder),'-dpng')

end
%% ========== Initial grain size ===================================

switch c.station
    case 'Miege'
        out = dlmread('./Input/Initial state/InitialGrainSizeProfile_Miege.csv');
    otherwise
        out = dlmread('./Input/Initial state/InitialGrainSizeProfile.csv');
end

depth = out(:,1); %old depth scale in m
old_graind = out(:,2); %density in kg/m^3 on the old real depth scale

[new_graind] = ConvertToGivenDepthScale(depth, old_graind, depth_act,'intensive');

graind_out = 2*ones(size(c.cdel));
graind_out(1:find(~isnan(new_graind),1,'last'),1,j) = new_graind(1:find(~isnan(new_graind),1,'last'));
   
if sum(isnan(graind_out))>0
        error('NaN in initial temperature values')
end
if c.verbose == 1
    f = figure('Visible','on','OuterPosition',[0 0 8 18]);
    scatter(depth,old_graind,'o')
    hold on
    stairs([0; depth_act],[graind_out; graind_out(end)], 'LineWidth',2)
    scatter([depth_act],[graind_out], 'LineWidth',2)
    legend('Original data','Model initial state','Location','NorthOutside')

    ylabel('Grain size (mm)')
    legend boxoff
    axis tight
    xlim([0 max(depth_act)])
    box on
    view([90 90])
        print(f,sprintf('%s/Initial_graind',c.OutputFolder),'-dpng')

end

%% ========== Initial water content ================================
switch c.station
    case 'Miege'
        [~, ind_12] = min(abs(depth_act-12));
        [~, ind_37] = min(abs(depth_act-37));
%         thick_aq = depth_act(ind_25)- depth_act(ind_15);

        pore_volume = c.cdel *1000 .* (1./rhofirn_ini - 1/c.rho_ice);
        pore_volume(snowc_ini<c.smallno) = 0;
        
        saturation = 1; %2000 / (sum(pore_volume(ind_12:ind_37))*c.rho_water);
        slwc_ini(ind_12:ind_37,:) = saturation * pore_volume(ind_12:ind_37);
        fprintf('Intital total lwc of %0.2f mm\n\n', ...
            1000*sum(slwc_ini(ind_12:ind_37,:)));
        
    case 'FA'
        filename = '.\RetMIP\Input files\lwc\RetMIP_lwc_FA.csv';
        delimiter = ';';
        startRow = 2;
        formatSpec = '%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        fclose(fileID);
        M = [dataArray{1:end-1}];
        clearvars filename delimiter startRow formatSpec fileID dataArray ans;
        
        depth = M(:,1);
        old_lwc = M(:,2);        
        [slwc_ini] = ConvertToGivenDepthScale(depth, old_lwc, depth_act,'intensive');
    otherwise
        slwc_ini = 0*depth_act;
        old_lwc = 0*depth;
end

if c.verbose == 1
    f = figure('Visible','on','OuterPosition',[0 0 8 18]);
    stairs(depth_act,slwc_ini.*1000,'r', 'LineWidth',2)
    if c.retmip == 1
        hold on
        scatter(depth,old_lwc.*1000, 'LineWidth',2)
    end
    
    hold on
    legend('Original data','Model initial state','Location','NorthOutside')
    ylabel('Initial water content (mm)')
    legend boxoff
    axis tight
    box on
%     title(sprintf('Initial grain size profile\n on %s',datestr(time_dt(1))))
    view([90 90])
    print(f,sprintf('%s/Initial_slwc.tif',c.OutputFolder),'-dtiff')
end

%% ========== add other initialization here ========================

end
