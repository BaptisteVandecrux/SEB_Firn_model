function [data_new] = ResampleTable(data,opt)
% ResampleData resamples a table every round hour.
%   data is the input table containing at least a field 'time' with Matlab
%   timestamps.
%
%   data_new is a table containing in 'time' the new time stamps (full 
%   hours) and all the other variables contained in data interpolated on
%   the new time time stamps.
%
if ~exist('opt','var')
    opt='linear';
end

 
    time_start = data.time(1);
    % finding full hour equal or preceding the first time stamp in
    % data.time
    if floor((time_start - floor(time_start))/(1/24)) ~= ...
        ((time_start - floor(time_start))/(1/24))

        time_start = floor(time_start) + floor((time_start - floor(time_start))/(1/24))*(1/24)+ 1/24;
    end

    time_end = data.time(end);
    % finding full hour equal or following the last time stamp in
    % data.time
    if floor((time_end - floor(time_end))/(1/24)) ~= ...
        ((time_end - floor(time_end))/(1/24))

        time_end = floor(time_end) + floor((time_end - floor(time_end))/(1/24))*(1/24) - 1/24;
    end

    % new time stamps
    new_time = [time_start:1/24:time_end]';
    VarNames = data.Properties.VariableNames;
    data_new = table;
    data_new.time = new_time;

    % interpolating all variables at the new time stamps
    for i = 1:size(data,2)
        if ~strcmp(VarNames{i},'time')
            if isnumeric(data.(VarNames{i}))
                interp_values = NaN(length(new_time), size(data.(VarNames{i}),2));
                for k = 1:size(data.(VarNames{i}),2)
                    interp_values(:,k) = interp1(data.time, data.(VarNames{i})(:,k), new_time,opt,'extrap');
                end
                data_new.(VarNames{i}) = interp_values;
                %if some data was interpolated between points of different origins,
                %we just round the interpolated origin to keep it integer
                if ~isempty(strfind(VarNames{i},'Origin'))
                    data_new.(VarNames{i}) = round(data_new.(VarNames{i}));
                end
            else
                for j = 1:length(new_time)
                    [~,ind] = min(abs(data.time-new_time(j)));
                    data_new.(VarNames{i})(j) = data.(VarNames{i})(ind);
                end
            end
        end
    end
end