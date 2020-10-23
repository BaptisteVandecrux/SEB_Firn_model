function [depth_act_from_surface, depth_act_fixed_origin, ...
    T_ice_mod, compaction, time_mod, snowbkt, rho_all,H_surf,...
    depth_thermistor, T_ice_obs, Surface_Height,tlwc] = ...
    LoadTempModel(filename,station)

%% Loading modelled temperature
    load(strcat(filename,'/run_param.mat'))

    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',filename,1);
    try Year = ncread(namefile,'Year');
    catch me
        namefile = sprintf('%s/%s_surf-bin-%i.nc',filename,station,1);
        Year = ncread(namefile,'Year');
    end
    Day = ncread(namefile,'Day');
    Hour = ncread(namefile,'Hour');
    snowbkt = ncread(namefile,'snowbkt');

    % extract subsurface variables
    namefile = sprintf('%s/subsurf-bin-%i.nc',filename,1);   
    if exist(namefile) == 2
        snowc = ncread(namefile,'snowc');
        snic = ncread(namefile,'snic');
    	slwc = ncread(namefile,'slwc');
        rho = ncread(namefile,'rho');
        T_ice_mod = ncread(namefile,'T_ice');
        namefile = sprintf('%s/subsurf-bin-%i.nc',filename,1);       
        compaction = ncread(namefile,'compaction');
    else
        snowc = ncread([filename '/snowc_bin_1.nc'],'snowc');
        snic = ncread([filename '/snic_bin_1.nc'],'snic');
        slwc = ncread([filename '/' station '_slwc_bin_1.nc'],'slwc');
        rho         = ncread([filename '/' station '_rho_bin_1.nc'],    'rho');
        T_ice_mod   = ncread([filename '/' station '_T_ice_bin_1.nc'],    'T_ice');
        compaction  = ncread([filename '/' station '_compaction_bin_1.nc'],'compaction');
    end

    tlwc = sum(slwc,1);
    % Update BV2017: (not used anymore) Liquid water does not participate to the actual thickness
    % of a layer. Ice does not participate as long as there is enough pore
    % space in the snow to accomodate it.
    pore_space = snowc .* c.rho_water.*( 1./rho - 1/c.rho_ice);
    excess_ice = max(0, snic * c.rho_water / c.rho_ice - pore_space);
    thickness_act = snowc.*(c.rho_water./rho) + excess_ice;

    % thickness_act = snowc.*(c.rho_water./rho) + snic.*(c.rho_water./c.rho_ice);

    depth_act=zeros(size(thickness_act));
    for j=1:length(thickness_act(1,:))
            depth_act(1:end,j)=cumsum(thickness_act(1:end,j));
    end

    rho_all= (snowc + snic)./...
                (snowc./rho + snic./c.rho_ice);
%     lwc = slwc(:,:) ./ thickness_act(:,:);

    time_mod = datenum(Year,1,Day,Hour,0,0);
%     TT = ones(c.jpgrnd,1) * time_mod';

    H_surf = depth_act(end,:)'-depth_act(end,1); %+snowbkt*1000/315;

    for i = 1:length(H_surf)-1
        if (H_surf(i+1)-H_surf(i))> c.new_bottom_lay-1
            H_surf(i+1:end) = H_surf(i+1:end) - c.new_bottom_lay*c.rho_water/c.rho_ice;
        end
    end    

    % depth scale
%     depth_act = vertcat(zeros(size(depth_act(1,:))), depth_act);
%     depth_act = vertcat(zeros(size(depth_act(1,:))), depth_act);
    depth_act_from_surface = depth_act;
    for i = 1:size(depth_act,1)
        depth_act(i,:) = depth_act(i,:) - H_surf' + H_surf(1);
    end
    depth_act_fixed_origin = depth_act;

    H_surf = depth_act_from_surface(end,:)'-depth_act_from_surface(end,1); %+snowbkt{i}*1000/315;

    for j = 1:length(H_surf)-1
        if (H_surf(j+1)-H_surf(j))> c.new_bottom_lay-1
            H_surf(j+1:end) = H_surf(j+1:end) - c.new_bottom_lay*c.rho_water/c.rho_ice;
        end
    end  
    
%% Exporting density
% Subsurface variables  
% data_subsurf = {rho_all};
% varname =  {'rho'};
% unit =  {'kg/m^3'};
% long_varname = {'Bulk firn density'};
% 
% for ii = 1:length(data_subsurf)
%     fprintf('Saving %s \n', varname{ii})
%     tic
%     WriteNC_2D(sprintf('%s/%s_%s_bin_1.nc',...
%         'C:\Data_save\2020 JoG data\Corrected\data out',...
%         c.station, varname{ii}), ...
%         time_mod, depth_act_from_surface, data_subsurf{ii}, 'Depth', 'm', ...
%         'Depth below the initial surface position',...
%         varname{ii}, unit{ii}, long_varname{ii});
%     toc
% end

%% Loading measured temperature
c.InputAWSFile = ['C:\Data_save\2020 JoG data\Corrected\AWS data\' c.InputAWSFile(26:end)];
        [~, ~, ~, ~, ~,...
    ~, ~, ~, ~, ~,~, ...
    ~, ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ~, ~,...
    ~, ~, ~, ~, T_ice_obs, ...
    depth_thermistor, Surface_Height, ~, ~, ~] = ...
    ExtractAWSData(c);
    T_ice_obs=T_ice_obs';
    depth_thermistor=depth_thermistor';

    
%     for i = 1:length(time_mod)
%         try [depth_thermistor(:,i), ind] = sort(depth_thermistor(:,i));
%         T_ice_obs(:,i) = T_ice_obs(ind,i);
%     end
end