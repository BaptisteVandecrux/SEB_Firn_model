function [perf_tot] = ModelPerformanceHandR(PerfFuncName, time_mod,time_obs,...
    H_surf, Surface_Height, depth_act, rho_all, fileID, c)
% Choice between:
% PerfFuncName = 'MSE';
% PerfFuncName = 'RMSE';
% PerfFuncName = 'SS';
% PerfFuncName = 'NashSut';
    
PerfFunc = str2func(PerfFuncName);
fprintf(fileID, '\n Performance of the model\n');

% ------- surface height -------
[perf_H_surf, n_H_surf] = PerfFunc(time_mod, H_surf-H_surf(1),time_mod, ...
    Surface_Height-Surface_Height(1));
fprintf(fileID, 'Surface height:\t %s = %0.3f \t n = %i\n',PerfFuncName,...
    perf_H_surf, n_H_surf);

% ------ Density Profiles ------
switch c.station
    case 'DYE-2'
        i_core = 12:14;
    case 'CP1'
        i_core = 15:27;
    case 'KAN-U'
        i_core = 1:11;
    otherwise
        NameStation = c.station;
        i_core = FindCore(Core,'NearestCodeLocation',NameStation);
end

 perf_dens = zeros(size(depth_obs,1),1);
 n_dens = zeros(size(depth_obs,1),1);
 for i = 1:length(i_core)
     
     depth = Core{i_core(i)}.Data.Depth;
     density = Core{i_core(i)}.Data.Density;
     A_obs = trapz (depth,density);
     
     [~, ind_core] = min(abs(time_mod,datenum(Core{i_core(i)}.Info.DateCored)));
     depth_mod = depth_act(:,ind_core);
     density_mod = rho_all(:,ind_core);
     A_mod = trapz(depth_mod,density_mod);
     
     perf_dens(i) = abs(A_mod - A_obs);

 end

 %Overall

perf_tot = sum( [perf_H_surf/0.5, perf_dens(:)'] .* w );

fprintf(fileID,'\n Standardized, weighted-averaged %s: %0.2f\n',PerfFuncName, perf_tot);

end