% Surface energy and mass budget model for ice sheets, by Dirk van As.
% The model can be run for a single location with a time series of p, T, RH, WS, SR, and LRin,
% or for a transect for which the input variables are height-dependent.
% The current version lacks:
% - reduction of sub-surface density by sub-surface melt
%
% This is model version 2009/12, with plotting bits taken out and minute column introduced and year-2000 and fixed T_ice and snowthick ini.
% 
% Update November 2015, by Robert S. Fausto (RSF)
% - HIRHAM5 subsurface scheme has replaced the former one. The subsurface scheme now includes retention by capilary forces and dry snow densification.
% - The subsurface subroutine is called "subsurface_hirham". See subroutine for description of parameters and variab<=s.
% - Layers mass (or water equivalent thickness) is constant in time. Each layer has a part which is snow (snowc), ice (snic) and water (slwc).
% - sum water equivalent thickness of layer n is thickness_weq(n) = snowc(n)+snic(n)+slwc(n). This number is constant in time.
% - Partitioning into snow, water and ice varies from time to time, and so does the density of the snow fraction (the variab<= rhofirn).
% - This means that the actual thickness (as opposed to water eqiv), <=tâ€™s call that â€?actâ€? (as opposed to â€œweqâ€?), is:
%   thickness_act(n) = snowc(n)*(rho_w/rhofirn(n)) + snic*(rho_w/c.rho_ice) + slwc
%
% Update Spring 2016 by Peter Langen, Robert Fausto
% - New percolation scheme translated from Peter Langen's work: possibility
% to choose between a standard bucket scheme (donoDarcy =1) and a bucket
% scheme that passes to the next layer only the amount that would be
% allowed by a Darcy flow (do_no_darcy = 0).
% - New superimposed ice formulation as translated from Peter Langen's
% FORTRAN version
%
% other updates 2016-2017 by Baptiste Vandecrux
% - constants and parameter defined in external file (see Input folder)
% - routines to set initial density/temperature profiles from data (see
% IniTemperatureDensity function)
% - Tdeep changed to the mean temperature of the elevation bin over the
% study period (see PrepareTransect function)
% - parameter C_k for the fraction of potential refreezing occuring (see
% refreeze function)
% - Choices between several parametrization for fresh snow density
% (including a WS dependant).
% - Choice between different definition of precipitation
% - Possibility to run it as conduction model for use as in Humphrey et al. 2012
% - Lefebvre et al. 2003 for the roughness scale of snow and ice
% - Calonne 2012 for permeability
% - Runoff according to a Darcy flow through saturated snow
% - variable layer thickness. Dynamic merging/splitting of layers.
% - initial layer thickness dependant on the accumulation

function [RunName, c] = HHsubsurf(param)
% disp('Running...')
tic

set(0,'defaultfigurepaperunits','centimeters');
   set(0,'DefaultAxesFontSize',15)
   set(0,'defaultfigurecolor','w');
set(0,'defaultfigureinverthardcopy','off');
set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);

% #### Constant definition ####
% All constant values are defined in a set of csv file in the Input folder.
% They can be modiefied there or by giving new values in the "param"
% variable. The values of the constant given in param will overright the
% ones extracted from the csv files. The fieldnames in param should be the
% same as is c.
c = ImportConst(param);
[RunName, c] = OutputName(c,c.station);
diary(sprintf('%s/log.txt',c.OutputFolder));

[time, year, day, hour, pres,...
    T1, T2, z_T1, z_T2, o_T1,o_T2, ...
    RH1, RH2, z_RH1, z_RH2, o_RH1, o_RH2, ...
    WS1, WS2, z_WS1, z_WS2, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_AWS, c] = ...
    ExtractAWSData(c);

[elev, pres, T, RH, WS, SRin, LRin, rho_atm, nu, ~, ~, Tdeep] ...
    = PrepareTransect(pres, T2, RH2, WS2, SRin, LRin, c);

[H_comp, GF, GFsurf, GFsubsurf, H_melt_weq, H_rain,...
H_subl, H_surf,H_surf_weq, H_snow, L, LHF, meltflux, meltflux_internal, ...
rho, ~, runoff, SHF, SRnet,SRout_mdl,  T_ice, grndc, grndd , ~, ~, ~, zrogl, ...
pdgrain, refreezing, snowbkt, theta_2m, q_2m, ws_10m, Tsurf, snowthick,...
z_icehorizon, Re,Ri, err]...
= IniVar(c);

%Initialization of freshsnow density for both precipitation at the surface
%and use in the subsurface scheme
c.rho_snow = IniRhoSnow(T1, WS1, elev, c);

% Converts RH from relative to water to relative to ice
% RH_wrt_w = RH;
% RH = RHwater2ice(RH_wrt_w, T, pres);

%Calculates specific humidity and saturation (needs RH relative to ice!)
[RH1, q1] = SpecHumSat(RH1, T1, pres, c);
[RH2, q2] = SpecHumSat(RH2, T2, pres, c);

% q = q*1000;
%Calculated precipitation types and temp
if sum(strcmp(data_AWS.Properties.VariableNames,'Snowfallmweq'))>0
    snowfall = data_AWS.Snowfallmweq;
    rainfall = zeros(size(T1));
    T_rain = zeros(size(T1));
else
    [snowfall, rainfall, T_rain, c] = ...
        Precipitation(time, T1, LRin_AWS,RH1, Surface_Height, c);
end

c = CalculateMeanAccumulation(time,snowfall, c);

% Update BV 2018
if c.track_density
    density_avg_20 = NaN(6,c.M);
    % figure
    % hold on
end
                        
%% START OF SPATIAL LOOP -----------------------------------------------------------------------
for j=1:c.elev_bins
    c.Tdeep = Tdeep(j);
% Initial temperature and density profile
%         if ~isnan(LRout(1))
%             Tsurf_ini = (LRout(1)/c.sigma).^(1/4);
%         else
        Tsurf(1,j) = T2(1);
%         end

    [T_ice, rhofirn,rho(:,1),snic, snowc, slwc, pdgrain(:,1)] = ...
        InitializationSubsurface(T_ice_obs, depth_thermistor, T_ice, ...
        time, Tsurf(1,j), j, c);
    
%     theta1 = T1(:,j) + z_T1(:,j) * c.g / c.c_pd;
    theta2 = T2(:,j) + z_T2(:,j) * c.g / c.c_pd;
    
    % START OF TIME LOOP -----------------------------------------------------------------------
    for k = 1:c.M

        %=========== Step 1/*: Update snowthickness and instrument heights ====================================
        [snowthick, z_icehorizon] = ...
            UpdateSnowThickness(snowthick,z_icehorizon, k, j, c);
   
        % ========== Step 3/*: shortwave radiation balance snow & ice penetration ====================================
        if k > 1
            rho(:,k) = rho(:,k-1);
        end
        [~, SRnet, T_ice, meltflux_internal, ~] = ...
            SRbalance (SRout, SRin, SRnet,...
                z_icehorizon, snowthick, T_ice, rho, j, k, c);

         % ========== Step 5/*:  Surface temperature calculation ====================================


        k_eff = 0.021 + 2.5e-6*rho(:,k).^2 ;        
        % effective conductivity by Anderson 1976, is ok at limits
        % thickness of the first layer in m weq for thermal transfer
        thick_first_lay = snic(1) + snowc(1);

        % starting surface temperature solving
            c.o_THF = 1;
            if and(c.solve_T_surf == 0, ~isnan(Tsurf_obs(k)))
                % if user asked for using measured surface temperature instead
                % of solving for it and that there is a measurement available
                iter_max_EB = 1; % we do one iteration in the solver bellow
                Tsurf(k,j) = Tsurf_obs(k); % we use measured surface temp
            else
                iter_max_EB = c.iter_max_EB; % otherwise just using standard value
            end
            
            % Prepare parameters needed for SEB
            EB_prev = 1;
            dTsurf = c.dTsurf_ini ;  % Initial surface temperature step in search for EB=0 (C)

% figure
% hold on

            for findbalance = 1 : iter_max_EB
                % SENSIBLE AND LATENT HEAT FLUX -------------------------------
                    err(k,j) = 0;

                        c.o_THF = 1;
                        [L(k,j), LHF(k,j), SHF(k,j), theta_2m(k,j), q_2m(k,j), ws_10m(k,j),Re(k,j)] ...
                            = SensLatFluxes_bulk (WS2(k,j), nu(k,j), q2(k,j), snowthick(k,j), ...
                            Tsurf(k,j), theta2(k,j), pres(k,j), rho_atm(k,j),  z_WS2(k,j), z_T2(k,j), z_RH2(k,j), ...
                            snowbkt, rho(1,k), c);

                % SURFACE ENERGY BUDGET ---------------------------------------

                        [meltflux(k,j), Tsurf(k,j), dTsurf, EB_prev, stop] ...
                            = SurfEnergyBudget (SRnet, LRin(k,j), Tsurf(k,j), k_eff,thick_first_lay, ...
                            T_ice(:,k,j), T_rain(k,j),...
                            dTsurf, EB_prev, SHF(k,j), LHF(k,j), rainfall(k,j),c);
% scatter(findbalance,Tsurf(k,j))
% xlim([0 iter_max_EB])
% title(c.o_THF)
% pause(0.001)
         

                if stop
                    break
                end
            end %end loop surface energy balance

            if iter_max_EB ~= 1 && ...
                (findbalance == c.iter_max_EB && abs(meltflux(k,j)) >= 10*c.EB_max)
                error('Problem closing energy budget')
            end
            clear findbalance
        
        % ========== Step 6/*:  Mass Budget ====================================
    [dH_melt_weq, H_melt_weq, dH_subl_weq, H_subl, H_snow, H_rain, snowthick] = ...
        MassBudget (meltflux, H_subl, H_melt_weq, H_snow, H_rain, LHF, ...
        snowfall, rainfall, snowthick, elev, j, k, c);

    % in the case of the conduction model, the mass budget is calculated as
    % follows
        if c.ConductionModel == 1
            smoothed_Surface_Height= smooth(Surface_Height,24*7);
            if k>1
                dSurface_Height= -(smoothed_Surface_Height(k) - smoothed_Surface_Height(k-1)); %in real m
            else
                dSurface_Height= 0;
            end
            if dSurface_Height<= 0
                % if the surface height increase, it means that snow is
                % falling
                dH_melt_weq = 0;
                snowfall(k,j) = -dSurface_Height*c.rho_snow(k,j)/c.rho_water; %in m weq
                dH_subl_weq = 0;
            else
                %else we just say it has sublimated (quick way to make the
                %matter disappear in the subsurface scheme)
                dH_melt_weq = 0; %in m weq
                dH_subl_weq = -dSurface_Height*rho(1, k)/c.rho_water;
                snowfall(k,j) = 0;
            end
            c.liqmax =0;
            c.calc_CLliq = 0;
            Tsurf(k,j) = ((LRout(k) - (1-c.em)*LRin(k)) /(c.em*c.sigma))^(1/4);
        end
        
        % ========== Step 7/*:  Sub-surface model ====================================
        GF(2:c.z_ice_max) = -k_eff(2:c.z_ice_max).*(T_ice(1:c.z_ice_max-1,k,j)-T_ice(2:c.z_ice_max,k,j))./c.dz_ice;
        GFsurf(k,j) =-(k_eff(1)) * (Tsurf(k,j)- T_ice(2,k,j)) / thick_first_lay;
%         grndhflx = GFsurf(k,j);       
        pTsurf = Tsurf(k,j);
        ptsoil_in = T_ice(:,k,j);
        zsn = snowfall(k,j) + dH_subl_weq;
        snmel = -dH_melt_weq;
        raind = rainfall(k,j);
        c.rho_fresh_snow = c.rho_snow(k,j);
        
        if k==1
            grndc =T_ice(:,k,j);
            grndd(:) =0;
        end

        [snowc, snic, slwc, ptsoil_out, zrfrz, rhofirn,...
            supimp, pdgrain, zrogl, ~, grndc, grndd, ~, grndhflx,...
            dH_comp, snowbkt, compaction, c] ...
            = subsurface(pTsurf, grndc, grndd, slwc, snic, snowc, rhofirn, ...
            ptsoil_in, pdgrain, zsn, raind, snmel, zrogl, Tdeep(j),...
            snowbkt,c);

        
                        % Update BV 2018
                        if c.track_density
                            density_avg_20(1,k) = c.rho_avg_aft_comp;
                            density_avg_20(2,k) = c.rho_avg_aft_snow;
                            density_avg_20(3,k) = c.rho_avg_aft_subl;
                            density_avg_20(4,k) = c.rho_avg_aft_melt;
                            density_avg_20(5,k) = c.rho_avg_aft_runoff;
                            density_avg_20(6,k) = c.rho_avg_aft_rfrz;
%                             if k>10 
%                                 plot([k-1 k-1+[1:6]./6], [density_avg_20(6,k-1) ; density_avg_20(:,k)],'--o')
%                                 pause
%                             end
                        end
                            

        T_ice(:,k,j)=ptsoil_out;
%         Tsurf(k,j)= ptsoil_out(1);
        % bulk density
        rho(:,k)= (snowc + snic)./...
            (snowc./rhofirn + snic./c.rho_ice);
        refreezing(:,k,j) = zrfrz + supimp;
        runoff(k,j) = zrogl;
        sublimation(k,j) = dH_subl_weq;
        z_icehorizon = floor(snowthick(k,j)/c.dz_ice);
        GFsubsurf(k,j) = grndhflx;
        snowbkt_out(k,j) = snowbkt;
        
        if k> 1
            H_surf_weq(k,j) = H_surf_weq(k-1,j) - (runoff(k,j)-runoff(k-1,j)) ...
                + snowfall(k,j) + rainfall(k,j) + dH_subl_weq;

            % Update BV2017: With the layer-conservative model, the surface height
            % can be calculated outside of the sub-surface scheme assuming that the
            % bottom of the collumn remains at constant depth

            % cumulative dry compaction
            H_comp(k,j) = H_comp(k-1,j) + dH_comp; %in real m


        end
        
        if(snowthick(k,j) < 0)
         snowthick(k,j)=0;
        end

        % MODEL RUN PROGRESS ----------------------------------------------
        if c.verbose == 1
        if (mod(k-1 , 24) == 0)
            fprintf('%.2f,day of the year: %i.\n',time(k), day(k)); % print daily (24) time progress for k being hourly
        end
        end

        %SAVING SOME SUBSURFACE VARIABLES --------------------------------------------
        if k==1
            sav. z_T = zeros(c.M,1);
            sav. slwc = zeros(c.jpgrnd,c.M);
            sav. snic = zeros(c.jpgrnd,c.M);
            sav. snowc = zeros(c.jpgrnd,c.M);
            sav. snowc = zeros(c.jpgrnd,c.M);
            sav. pdgrain = zeros(c.jpgrnd,c.M);
            sav. rhofirn = zeros(c.jpgrnd,c.M);
            sav. subsurf_compaction = zeros(c.jpgrnd,c.M);
            sav. o_THF = zeros(c.M,1);
        end
        
        sav. slwc(:,k) = slwc;
        sav. snic(:,k) = snic;
        sav. snowc(:,k) = snowc;
        sav. pdgrain(:,k) = pdgrain;
        sav. rhofirn(:,k) = rhofirn;
        sav. subsurf_compaction(:,k) = compaction;
        sav.z_T(k) = z_T2(k);
        sav.o_THF(k) = c.o_THF;
    end  % END OF TIME LOOP -----------------------------------------------------------------------
    
    rainHF = c.rho_water.*c.c_w(1).*rainfall.*c.dev./c.dt_obs.*(T_rain-Tsurf(:,j));
    
    %----------------------------------------------------------------------------------------------------
    % INTERPOLATION BACK TO ORIGINAL TIME RESOLUTION
    M = c.M/c.dev;
    if c.dev ~= 1
        ipfac2 = 1:(M+1);
        time        = interp1(time,ipfac2);
        pres(j,:)     = interp1(pres(j,:),ipfac2);
        T(j,:)      = interp1(T(j,:),ipfac2);
        RH(j,:)     = interp1(RH(j,:),ipfac2);
        WS(j,:)     = interp1(WS(j,:),ipfac2);
        SRin(j,:)     = interp1(SRin(j,:),ipfac2);
        SRout_mdl(j,:)  = interp1(SRout_mdl(j,:),ipfac2);
        LRin(j,:)     = interp1(LRin(j,:),ipfac2);
        Tsurf(j,:)    = interp1(Tsurf(j,:),ipfac2);
        L(j,:)      = interp1(L(j,:),ipfac2);
        GFsurf(j,:)   = interp1(GFsurf(j,:),ipfac2);
        SHF(j,:)      = interp1(SHF(j,:),ipfac2);
        LHF(j,:)      = interp1(LHF(j,:),ipfac2);
        rainHF(j,:)   = interp1(rainHF(j,:),ipfac2);
        meltflux(j,:)   = interp1(meltflux(j,:),ipfac2);
        meltflux_internal(j,:) =  interp1(meltflux_internal(j,:),ipfac2);
        %         H_aws(j,:)    = interp1(H_aws(j,:),ipfac2);
        H_surf(j,:)   = interp1(H_surf(j,:),ipfac2);
        H_melt_weq(j,:)   = interp1(H_melt_weq(j,:),ipfac2);
        H_snow(j,:)   = interp1(H_snow(j,:),ipfac2);
        H_subl(j,:)   = interp1(H_subl(j,:),ipfac2);
        runoff(j,:)   = interp1(runoff(j,:),ipfac2)*c.dev;
        rainfall(j,:)   = interp1(rainfall(j,:),ipfac2)*c.dev;
        snowfall(j,:)   = interp1(snowfall(j,:),ipfac2)*c.dev;
        snowthick(j,:)  = interp1(snowthick(j,:),ipfac2);
    end
    
    %% Writing data to net cdf

    % Surface variables
    namefilesurf = sprintf('%s/surf-bin-%i.nc',c.OutputFolder,j);
    
    data = {time,       year,       day,            hour,   SRout(:,j), ...
            c.em*c.sigma*Tsurf(:,j).^4-(1-c.em)*LRin(:,j), SHF(:,j), LHF(:,j), ...
            GFsurf(:,j),GFsubsurf(:,j),    rainHF(:,j),        meltflux(:,j),  ...
            H_surf(:,j),    H_surf_weq(:,j),    H_melt_weq(:,j),    ...
            H_subl(:,j),    H_comp(:,j),        runoff(:,j),    snowthick(:,j), ...
            snowfall(:,j),  rainfall(:,j),      SRout(:,j)./SRin(:,j), ...
            Tsurf(:,j)-c.T_0, L(:,j), sav.z_T, sublimation(:,j), snowbkt_out(:,j),...
            theta_2m(:,j), RHice2water( spechum2relhum(theta_2m(:,j), pres, q_2m(:,j),c),theta_2m, pres), ws_10m(:,j)};
    varname =  {'Time'  'Year'  'Day'   'Hour' 'SRout_mdl' ...
        'LRout_mdl' 'SHF' 'LHF' ...
        'GF' 'GFsubsurf' 'rainHF' 'meltflux'...
        'H_surf' 'H_surf_weq' 'H_melt' ...
        'H_subl' 'H_comp' 'runoff' 'snowthick'...
        'snowfall' 'rainfall' 'Albedo' 'Tsurf' 'L_Ob', 'H_T','sublimation','snowbkt',...
        'theta_2m','RH_2m_wrtw','ws_10m'};
    unit =  {'yr' 'yr' 'd' 'h' 'Wm-2' ...
        'Wm-2' 'Wm-2' 'Wm-2' ...
        'Wm-2' 'Wm-2' 'Wm-2' 'Wm-2' ...
        'm' 'm_weq' 'm_weq'...
        'm_weq' 'm' 'm_weq' 'm_weq'...
        'm_weq' 'm_weq' '-' 'C' 'm', 'm', 'm_weq','m_weq','K','%','ms-1'};
     for i=1:length(data)
        WriteNetCDF(namefilesurf,data{i},varname{i},1, c.M, unit{i});
     end
      clear data varname unit
      
    % Subsurface variables
    namefilesubsurf = sprintf('%s/subsurf-bin-%i.nc',c.OutputFolder,j);
    data = {T_ice(:,:,j) sav.rhofirn sav.slwc sav.snic sav.snowc sav.pdgrain...
        refreezing(:,:,j) sav.subsurf_compaction};
    varname =  {'T_ice' 'rho' 'slwc' 'snic' 'snowc' 'dgrain' 'rfrz' 'compaction'};
    unit =  {'C' 'kg/m^3' 'mweq' 'mweq' 'mweq' 'mweq' 'mm' 'm'};

    for i = 1:length(data)
        WriteNetCDF(namefilesubsurf, data{i}, varname{i}, c.M, c.jpgrnd, unit{i})
    end
    clear data varname unit
    
    if c.track_density
            % Subsurface variables
        namefilesubsurf = sprintf('%s/track_density.nc',c.OutputFolder);
        data = {density_avg_20};
        varname =  {'density_avg_20'};
        unit =  {'kg/m^3'};

        for i = 1:length(data)
            WriteNetCDF(namefilesubsurf, data{i}, varname{i}, c.M, 6, unit{i})
        end
        clear data varname unit
    end
  
end  % END OF SPATIAL LOOP -----------------------------------------------------------------------
save(strcat(c.OutputFolder,'/run_param.mat'),'c')

% Output the latent and sensible heat fluxes along with the stability
% metrics. Will be done in a smarter way in a future release.
% M = [time SHF LHF sav.o_THF Re Ri err];
% dlmwrite(sprintf('THF_%s_%i.csv',c.station ,c.THF_calc),M,'Delimiter',',','precision',9);

disp ('Done...')
toc

end



