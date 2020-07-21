function ts2=DailyAvgTS(ts)
    time =  floor(datenum(ts.TimeInfo.StartDate) + ts.Time);
    [unDates, ~, subs] = unique(time);

    % Accumulate by day
    data = accumarray(subs, ts.Data, [], @mean);
    time=unDates;
    ts2=timeseries(data,datestr(time));
end