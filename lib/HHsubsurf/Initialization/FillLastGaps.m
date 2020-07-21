function data_out = FillLastGaps(data_out,c)
    VarNames = {'AirPressurehPa' 'AirTemperatureC'  'AirTemperature2C' ...
        'RelativeHumidity' 'RelativeHumidity2' ...
        'WindSpeedms' 'WindSpeed2ms' ...
        'ShortwaveRadiationDownWm2' 'ShortwaveRadiationDown_CorWm2' ...
        'ShortwaveRadiationUpWm2' 'ShortwaveRadiationUp_CorWm2'...
        'LongwaveRadiationDownWm2' 'HeightWindSpeedm' 'HeightTemperaturem', ...
        'HeightHumiditym' , ...
         'HeightWindSpeed2m' 'HeightTemperature2m', ...
        'HeightHumidity2m'};

    for i =1:length(VarNames)
        if sum(strcmp(data_out.Properties.VariableNames,VarNames{i})) == 0
            ind_first_nonnan(i) = 1;
            ind_last_nonnan(i) = height(data_out);
            continue
        end       

        if isempty(strfind(VarNames{i},'Height'))
            [missingvalues, data_out.(VarNames{i})] = ...
                FillMissingValues(data_out.(VarNames{i}), 'linear');
            if c.verbose == 1 && sum(missingvalues)>1
                    fprintf('%s: %i missing values linearly interpolated\n',...
                        VarNames{i}, missingvalues);
            end
        else
            [missingvalues(i), data_out.(VarNames{i})] = ...
                FillMissingValues(data_out.(VarNames{i}), 'nearest','extrap');    
        end
        ind_first_nonnan(i) = ind_isnan(data_out.(VarNames{i}),'beginning');
        ind_last_nonnan(i) = ind_isnan(data_out.(VarNames{i}),'end');
    end

      % Removing dates for which data could not be interpolated (at the beginning
    % and at the end
    if max(ind_first_nonnan)>1 || min(ind_last_nonnan)<height(data_out)
        disp('Cropping file because of missing')
        disp('values at the begining or at the')
        disp('end of the record.')
        data_out = data_out(max(ind_first_nonnan):min(ind_last_nonnan),:);
    end
        disp('')
end