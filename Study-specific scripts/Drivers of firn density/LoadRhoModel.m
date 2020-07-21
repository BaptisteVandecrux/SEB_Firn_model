function [depth_act_from_surface, depth_act_fixed_origin, rho_all, time_mod, rho_20] = LoadRhoModel(filename)
    load(strcat(filename,'/run_param.mat'))

    % extract track density
    namefile = sprintf('%s/track_density.nc',filename);
    density_avg_20 = ncread(namefile,'density_avg_20');

    rho_20 = density_avg_20(1,:);
    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',filename,1);
    Year = ncread(namefile,'Year');
    Day = ncread(namefile,'Day');
    Hour = ncread(namefile,'Hour');

    % extract subsurface variables
    namefile = sprintf('%s/subsurf-bin-%i.nc',filename,1);       
    snowc = ncread(namefile,'snowc');
    snic = ncread(namefile,'snic');
    slwc = ncread(namefile,'slwc');
    rho = ncread(namefile,'rho');

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
    lwc = slwc(:,:) ./ thickness_act(:,:);

    time_mod = datenum(Year,1,Day,Hour,0,0);
    TT = ones(c.jpgrnd,1) * time_mod';

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
end