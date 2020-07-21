function [data_DJF] = AvgTableDJF(data,func)
    

if data.time(1)< 700000
    % then it is in decimal year
    out_dec_year =1 ;
    data.time = datenum(data.time,1,1);
else
    out_dec_year = 0;
end
    temp = datevec(data.time);
    ind_DJF = ismember(temp(:,2),[12 1 2]);
    data_DJF = data(ind_DJF,:);
 
    data_DJF = AvgTable(data_DJF,'water-yearly',func,90);
    temp = datevec(data_DJF.time);
    data_DJF.Year = temp(:,1);
   varname = data.Properties.VariableNames;
    % if the hourly data stops after the second of June then we
    % consider the summer not to be complete and we assign NaN
    if data.time(1)>datenum(temp(1,1),12,3)
        for i = 1:length(varname)
            if strcmp(varname{i},'time')
                continue
            else
                data_DJF.(varname{i})(1) = NaN;
            end
        end
    end
    % if the hourly data stops before the first of September then we
    % consider the summer not to be complete and we assign NaN
    if data.time(end)<datenum(temp(end,1),2,31)
        for i = 1:length(varname)
            if strcmp(varname{i},'time')
                continue
            else
                data_DJF.(varname{i})(end) = NaN;
            end
        end
    end
    
%     if out_dec_year==1
%         data_DJF.time = decyear(datestr(data_DJF.time),'dd-mmm-yyyy');
%     end
    
end