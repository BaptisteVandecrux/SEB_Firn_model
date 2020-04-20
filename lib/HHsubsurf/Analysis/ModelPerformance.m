function [perf_tot] = ModelPerformance(PerfFuncName, time_mod,time_obs, Tsurf, ...
    LRout, H_surf, depth_obs, Surface_Height, T_subsurf_mod, T_obs, ...
    fileID, c)
% Choice between:
% PerfFuncName = 'MSE';
% PerfFuncName = 'RMSE';
% PerfFuncName = 'SS';
% PerfFuncName = 'NashSut';
    
PerfFunc = str2func(PerfFuncName);
fprintf(fileID, '\n Performance of the model\n');
%----- surface temperature -----
if sum(~isnan(LRout))>10
    [perf_Tsurf, n_Tsurf] = PerfFunc(time_mod, Tsurf+c.T_0,time_mod, ...
        (LRout/c.sigma).^(1/4));
    fprintf(fileID, 'Surface temperature:\t\t\t %s = %0.3f \t n = %i\n',...
        PerfFuncName,perf_Tsurf, n_Tsurf);
else
    perf_Tsurf = 0;
    n_Tsurf =0;    
end
% ------- surface height -------
[perf_H_surf, n_H_surf] = PerfFunc(time_mod, H_surf-H_surf(1),time_mod, ...
    Surface_Height-Surface_Height(1));
fprintf(fileID, 'Surface height:\t %s = %0.3f \t n = %i\n',PerfFuncName,...
    perf_H_surf, n_H_surf);

% ------ Subsurf. temperature ------
if sum(sum(~isnan(T_subsurf_mod)))>10
     perf_T_subsurf = zeros(size(depth_obs,1),1);
     n_T_subsurf = zeros(size(depth_obs,1),1);
     for i = 1:size(depth_obs,1)
         [perf_T_subsurf(i), n_T_subsurf(i)] = PerfFunc(time_mod, ...
             T_subsurf_mod(i,:)', time_obs, T_obs(i,:)');
        fprintf(fileID, 'Subsurf. temperature (therm. %i): %s = %0.3f \t n = %i\n',...
            i, PerfFuncName, perf_T_subsurf(i), n_T_subsurf(i));
     end
else
     for i = 1:size(depth_obs,1)
         perf_T_subsurf(i)= 0;
         n_T_subsurf(i) =0;
     end
end
 %Overall
 w = [n_Tsurf/3, n_H_surf/3, n_T_subsurf(:)'/3/8]...
     / sum([n_Tsurf, n_H_surf, n_T_subsurf(:)']) ;
perf_tot = sum( [perf_Tsurf/1.1, perf_H_surf/0.5, perf_T_subsurf(:)'./0.05] .* w );

fprintf(fileID,'\n Standardized, weighted-averaged %s: %0.2f\n',PerfFuncName, perf_tot);

end