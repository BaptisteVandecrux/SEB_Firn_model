function [data_JJA] = AvgTableJJA(data,func)
    

if data.time(1)< 700000
    % then it is in decimal year
    out_dec_year =1 ;
    data.time = datenum(data.time,1,1);
else
    out_dec_year = 0;
end
    temp = datevec(data.time);
    ind_JJA = ismember(temp(:,2),[6 7 8]);
    data_JJA = data(ind_JJA,:);
 
    data_JJA = AvgTable(data_JJA,'yearly2',func,90);
    temp = datevec(data_JJA.time);
    data_JJA.Year = temp(:,1);
   varname = data.Properties.VariableNames;
    % if the hourly data stops after the second of June then we
    % consider the summer not to be complete and we assign NaN
    if data.time(1)>datenum(temp(1,1),6,3)
        for i = 1:length(varname)
            if strcmp(varname{i},'time')
                continue
            else
                data_JJA.(varname{i})(1) = NaN;
            end
        end
    end
    % if the hourly data stops before the first of September then we
    % consider the summer not to be complete and we assign NaN
    if data.time(end)<datenum(temp(end,1),8,31)
        for i = 1:length(varname)
            if strcmp(varname{i},'time')
                continue
            else
                data_JJA.(varname{i})(end) = NaN;
            end
        end
    end
    
%     if out_dec_year==1
%         data_JJA.time = decyear(datestr(data_JJA.time),'dd-mmm-yyyy');
%     end
    
end