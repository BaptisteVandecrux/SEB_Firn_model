function [data_out] = AvgTable(varargin)
% exemple: AvgTable(data, 'yearly','sum')

data = varargin{1};
period = varargin{2}; %period on which the average should be calculated
num_min_no_nan = 0; % minimum percentage of non-nan values needed for an average to be calculated
operation='sum';         % type of calculation needed, works with 'mean' and 'sum'
 
switch length(varargin)
    case 3
        operation = varargin{3};
    case 4
        operation = varargin{3};
        num_min_no_nan = varargin{4};
end
    time_start = datetime(datestr(data.time(1)));  
    time_end = datetime(datestr(data.time(end)));  
    time_new = datenum(time_start);
    i = 1;

switch period
    case 'monthly'
           while time_new(end) < datenum(time_end)
               time_new(i) = datenum(time_start.Year, time_start.Month + i-1, 1);
               i = i + 1;
           end
    case 'yearly'
        while time_new(end) + 360 < datenum(time_end)
           time_new(i) = datenum(time_start.Year + i-1, 1, 1);
           i = i + 1;
        end
    case 'yearly2'
        while time_new(end) + 360 < datenum(time_end)
           time_new(i) = datenum(time_start.Year + i-1, 1, 1);
           i = i + 1;
        end
    case 'seasonaly'
        temp = time_start.Month - [9, 12, 3, 6];
        temp(temp>0) = NaN;
        missing_until_next_season = -max(temp);
        time_new(i) = datenum(time_start.Year, ...
            time_start.Month + missing_until_next_season,1);
        i = i+1;

        %accumulation season from
        while time_new(end) + 350 < datenum(time_end)
            temp = datetime(datestr(time_new(i-1)));
           time_new(i) = datenum(temp.Year , temp.Month + 3, 1);
           i = i + 1;
        end
    case 'water-yearly'
        temp = time_start.Month - [9, 9+12];
        temp(temp>0) = NaN;
        missing_until_next_season = -max(temp);
        time_new(i) = datenum(time_start.Year, ...
            time_start.Month + missing_until_next_season,1);
        i = i+1;

        %accumulation season from
        while time_new(end) < datenum(time_end)
            temp = datetime(datestr(time_new(i-1)));
           time_new(i) = datenum(temp.Year , temp.Month + 12, 1);
           i = i + 1;
        end
    case 'Jun-yearly'
        temp = time_start.Month - [6, 6+12];
        temp(temp>0) = NaN;
        missing_until_next_season = -max(temp);
        time_new(i) = datenum(time_start.Year, ...
            time_start.Month + missing_until_next_season,1);
        i = i+1;

        %accumulation season from
        while time_new(end) < datenum(time_end)
            temp = datetime(datestr(time_new(i-1)));
           time_new(i) = datenum(temp.Year , temp.Month + 12, 1);
           i = i + 1;
        end
    case 'three-hourly'
       while time_new(end) < datenum(time_end)
           time_new(i) = datenum(time_start.Year, time_start.Month, ...
               time_start.Day+ time_start.Hour/24 + time_start.Minute/60/24 + (3/24) * (i-1));
           i = i + 1;
       end
    case 'hourly'
       while time_new(end) < datenum(time_end)
           time_new(i) = datenum(time_start.Year, time_start.Month, ...
               time_start.Day + time_start.Hour/24 + time_start.Minute/60/24 + (1/24) * (i-1));
           i = i + 1;
       end
    case 'daily'
       while time_new(end) < datenum(time_end)
           time_new(i) = datenum(time_start.Year, time_start.Month, floor(time_start.Day )+ i-1);
           i = i + 1;
       end
end

%    disp(datestr(time_new)); 
   data_out = array2table(repmat(time_new',1,size(data,2)));
   data_out.Properties.VariableNames = data.Properties.VariableNames;

    VarName = data.Properties.VariableNames;
for i = 1:size(data,2)
%     fprintf('%i/%i\n',i,size(data,2))
    if ~strcmp(VarName{i},'time')  
        ind_binned = discretize(data.time,time_new);
        ind_nan =isnan(ind_binned);
        temp = data.(VarName{i})(~ind_nan);
        ind_binned = discretize(data.time(~ind_nan),time_new);
        
        switch operation
            case 'sum'
            data_out.(VarName{i}) = accumarray(ind_binned,temp,...
                size(data_out.time),@sum);
            case 'nanmean'
            data_out.(VarName{i}) = accumarray(ind_binned,temp,...
                size(data_out.time),@nanmean);
            case 'mean'
            if num_min_no_nan == 0
                data_out.(VarName{i}) = accumarray(ind_binned,data.(VarName{i})(~ind_nan),...
                    size(data_out.time),@nanmean);
            else
                func_hdl = str2func(sprintf('nanmean_%i',num_min_no_nan));
                data_out.(VarName{i}) = accumarray(ind_binned,data.(VarName{i})(~ind_nan),...
                    size(data_out.time),func_hdl);
            end
        end
    end
end
    temp = datevec(data_out.time);
    temp = temp(:,1);

switch period
    case 'yearly'
            % if the hourly data stops after the 5th of Jan then we
        % consider the year not to be complete and we assign NaN
           varname = data_out.Properties.VariableNames;

        if data.time(1)>datenum(temp(1),1,5)
            for i = 1:length(varname)
                if strcmp(varname{i},'time')
                    continue
                else
                    data_out.(varname{i})(1) = NaN;
                end
            end
        end
        % if the hourly data stops before the 25th of December then we
        % consider the year not to be complete and we assign NaN
        if data.time(end)<datenum(temp(end),12,25)
            for i = 1:length(varname)
                if strcmp(varname{i},'time')
                    continue
                else
                    data_out.(varname{i})(end) = NaN;
                end
            end
        end

    otherwise
        data_out(end,:) = [];
end
end
    
        
            