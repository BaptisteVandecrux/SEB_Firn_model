function [time_out, data_out] = DailyAverage(time_in, data_in)
    time =  time_in-0.875;
    if time (1)-floor(time(1)) ~= 0
        data_in(time<floor(time(1)+1)) = [];
        time(time<floor(time(1)+1)) = [];
    end
    if time(end)-floor(time(end)) ~= 0
        data_in(time>floor(time(end))) = [];
        time(time>floor(time(end))) = [];
    end
    time = floor(time);
    [unDates, ~, subs] = unique(time);

    % Accumulate by day
    data_out = accumarray(subs, data_in, [], @mean);
    time_out=unDates+0.875;
end