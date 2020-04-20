% Surface energy and mass budget model for ice sheets, by Dirk van As.
% The model can be run for a single location with a time series of p, T, RH, WS, SR, and LRin,
% or for a transect for which the input variables are height-dependent.
% The current version lacks:
% - reduction of sub-surface density by sub-surface melt
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
    WS1, WS2, ~, z_WS2, o_WS1, o_WS2,...
    SRin, SRout, LRin, LRout, T_ice_obs, ...
    depth_thermistor, Surface_Height, Tsurf_obs, data_AWS, c] = ...
    ExtractAWSData(c);

if c.retmip
    [data_AWS,Tsurf_obs, pres, T1, T2, z_T1, z_T2, ...
        o_T1, o_T2,RH1, RH2, z_RH1, z_RH2, ...
        o_RH1, o_RH2,WS1, WS2, ~, z_WS2, ...
        o_WS1, o_WS2, SRin, SRout, LRin, LRout, c] = PrepareForRetMIP(c);
end

[elev, pres, ~, ~, ~, SRin, LRin, rho_atm, nu, ~, ~, Tdeep] ...
    = PrepareTransect(pres, T2, RH2, WS2, SRin, LRin, c);

%Initialization of freshsnow density for both precipitation at the surface
%and use in the subsurface scheme
c.rho_snow = IniRhoSnow(T1, WS1, elev, c);

% Converts RH from relative to water to relative to ice
% RH_wrt_w = RH;
% RH = RHwater2ice(RH_wrt_w, T, pres);

%Calculates specific humidity and saturation (needs RH relative to ice!)
[RH1, q1] = SpecHumSat(RH1, T1, pres, c);
[~, q2] = SpecHumSat(RH2, T2, pres, c);

[H_comp, GF, GFsurf, GFsubsurf, melt_mweq,snowbkt_out,...
sublimation_mweq, SMB_mweq,  ~, L, LHF, meltflux, ~, ...
rho,  runoff,...
SHF, SRnet,~,  T_ice, grndc, grndd , ...
pdgrain, refreezing,  theta_2m,q_2m,ws_10m, Tsurf, snowthick,...
z_icehorizon,Re,Ri,err]...
= IniVar(c);

%Calculated precipitation types and temp
if sum(strcmp(data_AWS.Properties.VariableNames,'Snowfallmweq'))>0
    snowfall = data_AWS.Snowfallmweq;
    if sum(strcmp(data_AWS.Properties.VariableNames,'Rainfallmweq'))>0
        rainfall=max(0,data_AWS.Rainfallmweq);
        T_rain = max(273.15,T2);
    else
        rainfall = zeros(size(T1));
        T_rain = zeros(size(T1));
    end
else
    [snowfall, rainfall, T_rain, c] = ...
        Precipitation(time, T1, LRin,RH1, Surface_Height, c);
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

    Tsurf(1,j) = mean(T2(1:24))-2.4;
% The 2 m air temperature and IR skin temperature are similar during peak 
% solar irradiance, with the mean difference in temperature equal to -0.32oC 
% when incoming solar radiation is greater than 600 W m2. There is a larger 
% difference between the two during the nighttime, with 2m air temperature
% higher than skin temperature by an average of 2.4??C when incoming 
% radiation is less than 200 Wm2. 
%  https://doi.org/10.5194/tc-12-907-2018

    [T_ice, rhofirn,rho(:,1),snic, snowc, slwc, pdgrain(:,1)] = ...
        InitializationSubsurface(T_ice_obs, depth_thermistor, T_ice, ...
        time, Tsurf(1,j), j, c);
    
    % preparing some variables
%     theta1 = T1(:,j) + z_T1(:,j) * c.g / c.c_pd;
    theta2 = T2(:,j) + z_T2(:,j) * c.g / c.c_pd;
%     theta_v1 = theta1 .* (1 + ((1 - c.es)/c.es).*q1);
    theta2_v = theta2 .* (1 + ((1 - c.es)/c.es).*q2);

    % Total precipitation in m of water
    if elev(j) <c.prec_cutoff
        snowfall(:,j) = snowfall(:,j)/2;
        rainfall(:,j) = rainfall(:,j)/2;
    end
    
    err(:,j) = 0;
    o_THF = err(:,j)+1;
    
    if c.THF_calc == 2 || c.THF_calc == 3
        % Testing conditions required for profile method
        % Weather variables from various origins
        test = o_T1(:,j) + o_T2(:,j) + o_RH1(:,j) + o_RH2(:,j) + o_WS1(:,j) + o_WS2(:,j) ~= 0;
        err(test,j) = 1;
            
        % WS below 1 m s-1
        test = or(WS1(:,j)<1, WS2(:,j)<1);
        err(test,j) = 2;   
            
        % Weather variables from same origin but T and RH measured at different height
        test = z_T1(:,j) ~= z_RH1(:,j);
        err(test,j) = 3;
            
        % WS2 <= WS1
        test = WS2(:,j)-WS1(:,j)<=0;
        err(test,j) = 4;
            
        if err(:,j) == 0       
            [LHF(:,j), SHF(:,j), theta_2m(:,j), q_2m(:,j), ws_10m(:,j),Ri(:,j)] ...
                = SensLatFluxes_profile (z_T1(:,j),z_T2(:,j),...
                T1(:,j),T2(:,j),q1(:,j),q2(:,j),WS1(:,j),WS2(:,j),pres(:,j), c);

            % Ri out of bound
            test =  isnan(Ri(:,j));
            err(test,j) = 5;
            
            % Unrealistically high SHF and LHF
            test =   abs(LHF(:,j))> 200 || abs(SHF(:,j))> 200;
            err(test,j) = 6;
            
            % Surface temperature below 0K
            test =   Tsurf(:,j)<=0 ;
            err(test,j) = 7;
            
            % Profile method output NaN (not from arleady listed errors)
            test =   isnan(LHF(:,j)) ||isnan(Tsurf(:,j)) || isnan(SHF(:,j)) ;
            err(test,j) = 8;
            
            LHF(err~=0,j)=NaN;
            SHF(err~=0,j)=NaN;
        end
        
        o_THF(err==0) = 2;
        o_THF(err~=0) = 1;
        
        if c.THF_calc == 3
            LHF2 = LHF;
            SHF2 = SHF;
            o_THF2 = o_THF;

            LHF(:)=NaN;
            SHF(:)=NaN;
            o_THF(:) = 1;
        end
    end
%     figure
    % START OF TIME LOOP -----------------------------------------------------------------------
    for k = 1:c.M
        %=========== Step 1/*: Update snowthickness and instrument heights ====================================
        [snowthick, z_icehorizon] = ...
            UpdateSnowThickness(snowthick,z_icehorizon, k, j, c);
        % ========== Step 3/*: shortwave radiation balance snow & ice penetration ====================================
        if k > 1
            rho(:,k) = rho(:,k-1);
            snowbkt_out(k,j) = snowbkt_out(k-1,j);
        end
        [~, SRnet, T_ice, ~ , ~] = ...
            SRbalance (SRout, SRin, SRnet,...
                z_icehorizon, snowthick, T_ice, rho, j, k, c);

         % ========== Step 5/*:  Surface temperature calculation ====================================

        k_eff = 0.021 + 2.5e-6*rho(:,k).^2 ;        
        % effective conductivity by Anderson 1976, is ok at limits
        % thickness of the first layer in m weq for thermal transfer
        thick_first_lay = snic(1) + snowc(1);
        
        % starting surface temperature solving
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

        if o_THF(k,j) == 1
            % Update BV2017: z_0 is calculated once outside the SEB loop
            if snowthick > c.smallno
                % if there is snow
                if snowbkt_out(k,j) > c.smallno
                    % fresh snow
                    z_0 = c.z0_fresh_snow;
                else
                    % old snow from Lefebre et al (2003) JGR
                    z_0 = max(c.z0_old_snow, ...
                        c.z0_old_snow + (c.z0_ice -c.z0_old_snow)*(rho(1,k) - 600)/(920 - 600));
                end
            else
                % ice roughness length
                z_0 = c.z0_ice;
            end
        end
    
        for findbalance = 1 : iter_max_EB
            % SENSIBLE AND LATENT HEAT FLUX -------------------------------
            if o_THF(k,j) == 1
                [L(k,j), LHF(k,j), SHF(k,j), theta_2m(k,j), q_2m(k,j), ws_10m(k,j),Re(k,j)] ...
                    = SensLatFluxes_bulk (WS2(k,j), nu(k,j), q2(k,j), snowthick(k,j), ...
                    Tsurf(k,j), theta2(k,j),theta2_v(k,j), pres(k,j), rho_atm(k,j),  ...
                    z_WS2(k,j), z_T2(k,j), z_RH2(k,j), z_0, c);
            end

            % SURFACE ENERGY BUDGET ---------------------------------------

            [meltflux(k,j), Tsurf(k,j), dTsurf, EB_prev, stop] ...
                = SurfEnergyBudget (SRnet, LRin(k,j), Tsurf(k,j), k_eff,thick_first_lay, ...
                T_ice(:,k,j), T_rain(k,j),...
                dTsurf, EB_prev, SHF(k,j), LHF(k,j), rainfall(k,j),c);
% scatter(findbalance,Tsurf(k,j))
% xlim([0 iter_max_EB])
% title(o_THF)
% pause(0.001)

                if iter_max_EB == 1
                    % if we are using surface temperature it might have been
                    % modified by SurfEnergyBudget. So we assign it again.
                    Tsurf(k,j) = Tsurf_obs(k); 
                end             

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
        % in mweq
        melt_mweq(k,j) = meltflux(k,j)*c.dt_obs/c.L_fus/c.rho_water;   
        sublimation_mweq(k,j) = LHF(k,j)*c.dt_obs/c.L_sub/c.rho_water; % in mweq
        % positive LHF -> deposition -> dH_subl positive


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
                melt_mweq(k,j) = 0;
                snowfall(k,j) = -dSurface_Height*c.rho_snow(k,j)/c.rho_water; %in m weq
                sublimation_mweq(k,j) = 0;
            else
                %else we just say it has sublimated (quick way to make the
                %matter disappear in the subsurface scheme)
                melt_mweq(k,j) = 0; %in m weq
                sublimation_mweq(k,j) = -dSurface_Height*rho(1, k)/c.rho_water;
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
        zsn = snowfall(k,j) + sublimation_mweq(k,j);
        snmel = melt_mweq(k,j);
        raind = rainfall(k,j);
        c.rho_fresh_snow = c.rho_snow(k,j);
        
        if c.retmip
            zsn = data_AWS.acc_subl_mmweq(k)/1000;
            snmel = data_AWS.melt_mmweq(k)/1000;
        end
        
        if k==1
            grndc =T_ice(:,k,j);
            grndd(:) =0;
        end
        if strcmp(c.station,'Miege')
            [slwc] = MimicAquiferFlow(snowc, rhofirn, snic, slwc, k,  c);
        end
        
        [snowc, snic, slwc, T_ice(:,k,j), zrfrz, rhofirn,...
            supimp, pdgrain, runoff(k,j), ~, grndc, grndd, ~, GFsubsurf(k,j),...
            dH_comp, snowbkt_out(k,j), compaction, c] ...
            = subsurface(pTsurf, grndc, grndd, slwc, snic, snowc, rhofirn, ...
            ptsoil_in, pdgrain, zsn, raind, snmel,  Tdeep(j),...
            snowbkt_out(k,j),c);

        % Update BV 2018
        if c.track_density
            density_avg_20(1,k) = c.rhoCC20_aft_comp(1);
            density_avg_20(2,k) = c.rhoCC20_aft_snow(1);
            density_avg_20(3,k) = c.rhoCC20_aft_subl(1);
            density_avg_20(4,k) = c.rhoCC20_aft_melt(1);
            density_avg_20(5,k) = c.rhoCC20_aft_runoff(1);
            density_avg_20(6,k) = c.rhoCC20_aft_rfrz(1);
            CC20(1,k) = c.rhoCC20_aft_comp(2);
            CC20(2,k) = c.rhoCC20_aft_snow(2);
            CC20(3,k) = c.rhoCC20_aft_subl(2);
            CC20(4,k) = c.rhoCC20_aft_melt(2);
            CC20(5,k) = c.rhoCC20_aft_runoff(2);
            CC20(6,k) = c.rhoCC20_aft_rfrz(2);
        end
        
        % bulk density
        rho(:,k)= (snowc + snic)./...
            (snowc./rhofirn + snic./c.rho_ice);
        refreezing(:,k,j) = zrfrz + supimp;
        z_icehorizon = floor(snowthick(k,j)/c.dz_ice);
        
        if k> 1
            SMB_mweq(k,j) =  snowfall(k,j) - runoff(k,j) ...
                + rainfall(k,j) + sublimation_mweq(k,j);

            % Update BV2017: With the layer-conservative model, the surface height
            % can be calculated outside of the sub-surface scheme assuming that the
            % bottom of the collumn remains at constant depth

            % cumulative dry compaction
            H_comp(k,j) = H_comp(k-1,j) + dH_comp; %in real m
        end
        
        if(snowthick(k,j) < 0)
         snowthick(k,j)=0;
        end

        % for the conduction model the temperature profile can be resetted
        % at fixed interval
        if c.ConductionModel == 1
            if (mod(k-1, 24) == 0)
                if sum(~isnan(T_ice_obs(k,:)))>0
                    [Tsurf(k,j), T_reset] = ...
                        ResetTemp(depth_thermistor, LRin, LRout, T_ice_obs, ...
                        rho, T_ice,time, k, c);
%                     figure
                    depth_act = cumsum(c.cdel .*c.rho_water ./rho(:,k));
                    depth_act = [0; depth_act];

%                     scatter(depth_thermistor(k,depth_thermistor(k,:)~=0),...
%                         T_ice_obs(k,depth_thermistor(k,:) ~= 0), 'o')
%                     hold on
%                     stairs(depth_act(1:end-1),T_ice(:,k,j)-c.T_0)
                    
                T_ice(~isnan(T_reset),k,j) = T_reset(~isnan(T_reset));
%                     stairs(depth_act(1:end-1),T_ice(:,k,j)-c.T_0)
%                     legend('data','before reset','after reset','Location','South')
%                     xlabel('Depth (m)')
%                     ylabel('Temperature (deg C)')
%                     title(sprintf('%s',datestr(datenum(time(k),0,0))))
%                     view([90 90])
   
                    [zso_capa, zso_cond] = ice_heats (c);
                    [grndc, grndd, ~, ~]...
                        = update_tempdiff_params (rho(:,k), Tdeep(j)                    ...
                        , snowc, snic, T_ice(:,k,j), zso_cond, zso_capa, c);
                end
            end            
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
        end
        
        sav. slwc(:,k) = slwc;
        sav. snic(:,k) = snic;
        sav. snowc(:,k) = snowc;
        sav. pdgrain(:,k) = pdgrain;
        sav. rhofirn(:,k) = rhofirn;
        sav. subsurf_compaction(:,k) = compaction;
        sav.z_T(k) = z_T2(k);
    end  % END OF TIME LOOP -----------------------------------------------------------------------
    
    rainHF = c.rho_water.*c.c_w(1).*rainfall./c.dt_obs.*(T_rain-Tsurf(:,j));
    

%% Processing few variables
    thickness_act = sav.snowc.*(c.rho_water./sav.rhofirn )+ ...
        sav.snic .*(c.rho_water/c.rho_ice);
    depth_act = cumsum(thickness_act, 1);
    
    H_surf = depth_act(end,:)'-depth_act(end,1)+snowbkt_out(k,j)*1000/315;
    for i = 1:length(H_surf)-1
        if (H_surf(i+1)-H_surf(i))> c.new_bottom_lay-1
            H_surf(i+1:end) = H_surf(i+1:end) - c.new_bottom_lay*c.rho_water/c.rho_ice;
        end
    end
    if c.retmip
        meltflux(:,j) = data_AWS.melt_mmweq/1000*c.L_fus*c.rho_water/c.dt_obs; 
        snowfall(:,j) = max(0,data_AWS.acc_subl_mmweq/1000); 
        sublimation_mweq(:,j) = min(0,data_AWS.acc_subl_mmweq/1000); 
    end
    
    %% Writing data to net cdf
    data_surf = {year,       day,            hour,   LRin(:,j), ...
            c.em*c.sigma*Tsurf(:,j).^4-(1-c.em)*LRin(:,j), SHF(:,j), LHF(:,j), ...
            GFsurf(:,j),    rainHF(:,j),        meltflux(:,j),  ...
            H_surf(:,j),    SMB_mweq(:,j),    melt_mweq(:,j),    ...
            sublimation_mweq(:,j),    H_comp(:,j),        runoff(:,j),    snowthick(:,j), ...
            snowfall(:,j),  rainfall(:,j),      SRin(:,j), SRout(:,j), ...
            Tsurf(:,j), sav.z_T, snowbkt_out(:,j),...
            theta_2m(:,j), RHice2water( spechum2relhum(theta_2m(:,j),...
            pres, q_2m(:,j),c),theta_2m, pres), ws_10m(:,j)};

    data_subsurf = {T_ice(:,:,j) sav.rhofirn sav.slwc sav.snic sav.snowc sav.pdgrain...
        refreezing(:,:,j) sav.subsurf_compaction};

try WritingModelOutput(time,data_surf,depth_act, data_subsurf,j,  c)
catch me 
ajf = 0 ;
end
end  % END OF SPATIAL LOOP -----------------------------------------------------------------------

save(strcat(c.OutputFolder,'/run_param.mat'),'c')

if c.THF_calc == 3
    M = [time SHF LHF o_THF SHF2 LHF2 o_THF2 Re Ri err];
    dlmwrite(sprintf('./Output/THF study/THF_%s_%i.csv',c.station ,c.THF_calc),M,'Delimiter',',','precision',9);
    % disp ('Done...')
end
    toc

end



