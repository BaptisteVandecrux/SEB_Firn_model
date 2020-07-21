function [data_SON] = AvgTableSON(data,func)
    

if data.time(1)< 700000
    % then it is in decimal year
    out_dec_year =1 ;
    data.time = datenum(data.time,1,1);
else
    out_dec_year = 0;
end
    temp = datevec(data.time);
    ind_SON = ismember(temp(:,2),[9 10 11]);
    data_SON = data(ind_SON,:);
 
    data_SON = AvgTable(data_SON,'yearly2',func,90);
    temp = datevec(data_SON.time);
    data_SON.Year = temp(:,1);
   varname = data.Properties.VariableNames;
    % if the hourly data stops after the second of June then we
    % consider the summer not to be complete and we assign NaN
    if data.time(1)>datenum(temp(1,1),9,3)
        for i = 1:length(varname)
            if strcmp(varname{i},'time')
                continue
            else
                data_SON.(varname{i})(1) = NaN;
            end
        end
    end
    % if the hourly data stops before the first of September then we
    % consider the summer not to be complete and we assign NaN
    if data.time(end)<datenum(temp(end,1),11,31)
        for i = 1:length(varname)
            if strcmp(varname{i},'time')
                continue
            else
                data_SON.(varname{i})(end) = NaN;
            end
        end
    end
    
%     if out_dec_year==1
%         data_SON.time = decyear(datestr(data_SON.time),'dd-mmm-yyyy');
%     end
    
end