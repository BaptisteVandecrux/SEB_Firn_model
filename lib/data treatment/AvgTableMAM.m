function [data_MAM] = AvgTableMAM(data,func)
    

if data.time(1)< 700000
    % then it is in decimal year
    out_dec_year =1 ;
    data.time = datenum(data.time,1,1);
else
    out_dec_year = 0;
end
    temp = datevec(data.time);
    ind_MAM = ismember(temp(:,2),[3 4 5]);
    data_MAM_in = data(ind_MAM,:);
 
    data_MAM = AvgTable(data_MAM_in,'yearly2',func);
    temp = datevec(data_MAM.time);
    data_MAM.Year = temp(:,1);
   varname = data.Properties.VariableNames;
    % if the hourly data stops after the second of March then we
    % consider the summer not to be complete and we assign NaN
    if data.time(1)>datenum(temp(1,1),3,3)
        for i = 1:length(varname)
            if strcmp(varname{i},'time')
                continue
            else
                data_MAM.(varname{i})(1) = NaN;
            end
        end
    end
    % if the hourly data stops before the first of May then we
    % consider the summer not to be complete and we assign NaN
    if data.time(end)<datenum(temp(end,1),5,31)
        for i = 1:length(varname)
            if strcmp(varname{i},'time')
                continue
            else
                data_MAM.(varname{i})(end) = NaN;
            end
        end
    end
    
%     if out_dec_year==1
%         data_MAM.time = decyear(datestr(data_MAM.time),'dd-mmm-yyyy');
%     end
    
end