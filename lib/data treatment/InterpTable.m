function [data_out] = InterpTable(data, gap_max_length)
% interpolate a table filling only
VarList = {'AirTemperature2mC', 'AirTemperature1C', 'AirTemperature2C',...
    'ShortwaveRadiationDownWm2', 'ShortwaveRadiationUpWm2', 'NetRadiationWm2', ...
    'AirPressurehPa', 'RelativeHumidity1Perc',...
    'RelativeHumidity2Perc', 'WindSpeed1ms', 'WindSpeed2ms', ...
    'WindDirection1deg', 'WindDirection2deg','Data'};
%    disp(datestr(time_new)); 
    data_out = data;
    VarName = data.Properties.VariableNames;
    for i = 1:size(data,2)
%         fprintf('%i/%i\n',i,size(data,2))
        if sum(strcmp(VarList,VarName{i}))>0 && sum(~isnan(data.(VarName{i}))) > 10
           data_out.(VarName{i}) = interp1gap(data.(VarName{i}) ,...
                    gap_max_length,'linear')';

        end
        
    end
end