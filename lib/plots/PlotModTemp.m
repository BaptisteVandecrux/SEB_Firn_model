list ={'KAN-U_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10_7', ...
    'KAN-U_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10_2', ...
    'KAN-U_0_SiCL_pr0.001_Ck1.00_darcy_wh0.10_9'};
   vis = 'on';
f = figure('Visible', vis,'Outerposition',[1 1 20 15]);
ha = tight_subplot(3,1,0.05,[0.15 0.02],0.1);    
for kk = 1:3
    RunName = list{kk};
    OutputFolder = sprintf('./Output/%s',RunName);

    % extract run parameters
    load(strcat(OutputFolder,'/run_param.mat'))
    c.OutputFolder = OutputFolder;
    % extract surface variables
    namefile = sprintf('%s/surf-bin-%i.nc',OutputFolder,1);
    finfo = ncinfo(namefile);
    names={finfo.Variables.Name};
    for i= 1:size(finfo.Variables,2)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end

    % extract subsurface variables
    namefile = sprintf('%s/subsurf-bin-%i.nc',OutputFolder,1);
    finfo = ncinfo(namefile);
    names={finfo.Variables.Name};
    for i= 1:size(finfo.Variables,2)
        eval(sprintf('%s = ncread(''%s'',''%s'');', char(names{i}), namefile,char(names{i})));
    end
    fprintf('\nData extracted from nc files.\n');


    thickness_weq = snowc + snic +slwc;

    % Update BV2017: (not used anymore) Liquid water does not participate to the actual thickness
    % of a layer. Ice does not participate as long as there is enough pore
    % space in the snow to accomodate it.
    pore_space = snowc .* c.rho_water.*( 1./rho - 1/c.rho_ice);
    excess_ice = max(0, snic * c.rho_water / c.rho_ice - pore_space);
    thickness_act = snowc.*(c.rho_water./rho) + excess_ice;

    % thickness_act = snowc.*(c.rho_water./rho) + snic.*(c.rho_water./c.rho_ice);

    depth_act=zeros(size(thickness_act));
    depth_weq=zeros(size(thickness_weq));
    for j=1:length(thickness_act(1,:))
            depth_act(1:end,j)=cumsum(thickness_act(1:end,j));
            depth_weq(1:end,j)=cumsum(thickness_weq(1:end,j));
    end

    rho_all= (snowc + snic)./...
                (snowc./rho + snic./c.rho_ice);
    lwc = slwc(:,:) ./ thickness_act(:,:);

    time_mod = datenum(Year,1,Day,Hour,0,0);
    TT = ones(c.jpgrnd,1) * time_mod';

    H_surf_old = H_surf;
    H_surf = depth_act(end,:)'-depth_act(end,1); %+snowbkt*1000/315;

    for i = 1:length(H_surf)-1
        if (H_surf(i+1)-H_surf(i))> c.new_bottom_lay-1
            H_surf(i+1:end) = H_surf(i+1:end) - c.new_bottom_lay*c.rho_water/c.rho_ice;
        end
    end    

    if ~isfield(c,'verbose')
        c.verbose = 1;
    end
    % extracting observed

    [time_yr, year, day, hour, pres,...
        T1, T, z_T1, H_instr_temp, o_T1,o_T2, ...
        RH1, RH, z_RH1, H_instr_hum, o_RH1, o_RH2, ...
        WS1, WS, z_WS1, H_instr_wind, o_WS1, o_WS2,...
        SRin, SRout, LRin, LRout, T_ice_obs, ...
        depth_thermistor, Surface_Height, Tsurf_obs, data_AWS, c] = ...
        ExtractAWSData(c);

    time_obs = datenum(year,1,day,hour,0,0);
    depth_obs = depth_thermistor';
    depth_obs(depth_obs==0) = NaN;

    T = T' + c.T_0;
    T(isnan(depth_obs)) = NaN;
    TT_obs= repmat(time_obs',size(depth_obs,1),1);
    % depth scale
    depth_act_save = depth_act;
    depth_act = vertcat(zeros(size(depth_act(1,:))), depth_act);
    depth_weq = vertcat(zeros(size(depth_weq(1,:))), depth_weq);

    for i = 1:size(depth_act,1)
        depth_act(i,:) = depth_act(i,:) - H_surf' + H_surf(1);
    end

    %% ============ Modelled Temperature =========================
    disp('Plotting modelled temperature ')
    range_temp = -30:1:0;
    set(0,'defaultfigurecolor','w');


    set(f,'CurrentAxes',ha(kk))
    col = PlotTemp(TT,depth_act,T_ice - c.T_0,...
        'Interp', 'on',...
        'PlotTherm', 'no',...
        'PlotIsoTherm', 'no',...
        'DataIsoTherm',rho_all,...
        'ValueIsoTherm', [800 800],...
        'ShowLegend','no',...
        'XLabel','Date',...
        'YLabel','Depth (m)',...
        'CLabel','Modelled Temperature (^oC)',...
        'cmap','jet',...
        'Range', range_temp);
    plot(time_mod,-H_surf+H_surf(1),'LineWidth',2)
    ylim([min([-3; -H_surf+H_surf(1)]) 30])
    col.TickLength=col.TickLength*5;
    col.Position = [ .88    0.3300    0.0430    0.6000];
    temp = col.YTickLabel;
    for i = 1:length(temp)
        if i/2==floor(i/2)
            temp(i,:)=' ';
        end
    end
    col.YTickLabel = temp;
    col.Units = 'normalized';
        plot(time_mod([1 end]), [2.2 2.2],'k','LineWidth',1.5)

    if kk==3

    col.Position = [0.85    ha(3).Position(2)    0.0430    ha(1).Position(3)+ha(1).Position(4)-0.15];
    else
    col.Position = [2    ha(3).Position(2)    0.0430    ha(1).Position(3)+ha(1).Position(4)];
    end
    ylim([0 6])
    set(gca,'Color',[0.95 0.95 0.95]);
    % set(gca,'XMinorTick','on','YMinorTick','on')
    xlim([time_mod(24*25) time_mod(24*30*3)])
    set_monthly_tick(time_mod,gca);
    datetick('x','mmm-yy','keepticks','keeplimits')
    set(gca,'XTickLabelRotation',0)
    ha(kk).LineWidth = 1.1;


end
ha(1).XTickLabel = '';
ha(2).XTickLabel = '';
ylabel(ha(1),'')
ylabel(ha(3),'')
xlabel(ha(1),'')
xlabel(ha(2),'')
annotation(f,'textbox',...
        [0.11 0.77 0.04 0.05],...
        'String','a)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','off');
annotation(f,'textbox',...
        [0.11 0.5 0.04 0.05],...
        'String','b)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','off');
annotation(f,'textbox',...
        [0.11 0.23 0.04 0.05],...
        'String','c)',...
        'LineStyle','none',...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','off');
    colormap jet
print(f, sprintf('%s/KANU_hetero',OutputFolder), '-dpng','-r300')
if strcmp(vis,'off')
    close(f);
end
