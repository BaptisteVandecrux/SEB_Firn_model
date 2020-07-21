function [SS, n] = SS(time_mod, modelled, time_obs, observed)

if isdatetime(time_mod)
    time_mod = datenum(time_mod);
end
if isdatetime(time_obs)
    time_obs = datenum(time_obs);
end

%first need to resample and interpolate to values at same time steps
ts1 = timeseries(modelled, time_mod);
ts2 = timeseries(observed, time_obs);

[ts1 ts2] = synchronize(ts1,ts2,'union');

modelled = ts1.Data;
observed = ts2.Data;

%only comparing data available in both vectors
ind = and(~isnan(modelled) , ~isnan(observed));

% Sum of square
SS =sum( (modelled(ind) - observed(ind)).^2 );
% we remove the number of interpolated points when counting element of
% comparison
n = sum(ind) - (length(ts1.Time) - max(length(time_obs),length(time_mod)));
if isnan(SS)
    SS = 0;
    n =0;
end

end