function [data_out] = AvgTable(varargin)
% exemple: AvgTable(data, 'yearly','sum')

data = varargin{1};
period = varargin{2}; %period on which the average should be calculated
num_min_no_nan = 0; % minimum percentage of non-nan values needed for an average to be calculated
fcn=@sum;         % type of calculation needed, works with 'mean' and 'sum'
 
switch length(varargin)
    case 3
        fcn = varargin{3};
    case 4
        fcn = varargin{3};
        num_min_no_nan = varargin{4};
end
    time_start = datevec(data.time(1));  
    time_end = datevec(data.time(end));  
    time_new = data.time(1);
    i = 1;

switch period
    case 'monthly'
           while time_new(end) < datenum(time_end)
               time_new(i) = datenum(time_start(:,1), time_start(:,2) + i-1, 1);
               i = i + 1;
           end
    case 'yearly'
        while time_new(end) + 360 < datenum(time_end)
           time_new(i) = datenum(time_start(:,1) + i-1, 1, 1);
           i = i + 1;
        end
    case 'yearly2'
        while time_new(end) + 360 < datenum(time_end)+6*31
           time_new(i) = datenum(time_start(:,1) + i-1, 1, 1);
           i = i + 1;
        end
    case 'seasonaly'
        temp = time_start(:,2) - [9, 12, 3, 6];
        temp(temp>0) = NaN;
        missing_until_next_season = -max(temp);
        time_new(i) = datenum(time_start(:,1), ...
            time_start(:,2) + missing_until_next_season,1);
        i = i+1;

        %accumulation season from
        while time_new(end) + 350 < datenum(time_end)
            temp = datevec(time_new(i-1));
           time_new(i) = datenum(temp(:,1) , temp(:,2) + 3, 1);
           i = i + 1;
        end
    case 'water-yearly'
        temp = time_start(:,2) - [9, 9+12];
        temp(temp>0) = NaN;
        missing_until_next_season = -max(temp);
        time_new(i) = datenum(time_start(:,1), ...
            time_start(:,2) + missing_until_next_season,1);
        i = i+1;

        %accumulation season from
        while time_new(end) < datenum(time_end)
            temp = datevec(time_new(i-1));
           time_new(i) = datenum(temp(:,1) , temp(:,2) + 12, 1);
           i = i + 1;
        end
    case 'Jun-yearly'
        temp = time_start(:,2) - [6, 6+12];
        temp(temp>0) = NaN;
        missing_until_next_season = -max(temp);
        time_new(i) = datenum(time_start(:,1), ...
            time_start(:,2) + missing_until_next_season,1);
        i = i+1;

        %accumulation season from
        while time_new(end) < datenum(time_end)
            temp = datevec(time_new(i-1));
           time_new(i) = datenum(temp(:,1) , temp(:,2) + 12, 1);
           i = i + 1;
        end
    case 'three-hourly'
       while time_new(end) < datenum(time_end)
           time_new(i) = datenum(time_start(:,1), time_start(:,2), ...
               time_start(:,3)+ time_start(:,4)/24 + time_start(:,5)/60/24 + (3/24) * (i-1));
           i = i + 1;
       end
    case 'hourly'
       while time_new(end) < datenum(time_end)
           time_new(i) = datenum(time_start(:,1), time_start(:,2), ...
               time_start(:,3) + time_start(:,4)/24 + time_start(:,5)/60/24 + (1/24) * (i-1));
           i = i + 1;
       end
    case 'daily'
       while time_new(end) < datenum(time_end)
           time_new(i) = datenum(time_start(:,1), time_start(:,2), floor(time_start(:,3) )+ i-1);
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
        
        data_out.(VarName{i}) = accumarray(ind_binned,temp,...
            size(data_out.time),fcn);

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
    
        
            