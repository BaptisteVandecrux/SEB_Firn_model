function [] = WriteNC_1D(namefile, time, data, varname, unit, long_varname)
% This is an attempt to create a robust netcdf writting script.
% It adds increment to the prescribed file name when file already exists or
% is not compatible with the variable that is being written (f.e. when
% dimensions do not match with existing file).
% The scripts can add Unit, LongVariableName and Comment as variable
% attribute.
% Example:
% data = {sin([1:1000]./6) exp(1:1000) -1000:-1};
% varname = {'var1','var2','var3'};
% long_varname = {'long_var1','long_var2','long_var3'};
% unit = {'m' 's' '-'};
% time =datenum(2019,1,1,1:1000,0,0);
% 
% WriteNC_1D('example.nc', time, data, varname, unit, long_varname)
% Baptiste Vandecrux bav@geus.dk
% ====================================================================

if time(1)<3000
    % then it is decimal year
    time = datenum(time,1,1) - datenum(1900,1,1,0,0,0);
elseif time(1)>401767
    % then it is in Matlab format
    time = time - datenum(1900,1,1,0,0,0);
else
    % then it is already in the good format
    warning('Time vector already in days since 1900-1-1 0:0:0')
end

if isrow(time)
    time = time';
end
     WriteVarToNC(namefile,time, ...
         'VariableName','time', ...
         'DimensionNames',{'time'},...
         'Unit', 'days since 1900-1-1 0:0:0', ...
         'LongVariableName','Time (UTC)');


         for i=1:length(data)
                if isrow(data{i})
                    data{i} = data{i}';
                end
                 WriteVarToNC(namefile,data{i}, ...
                     'VariableName',varname{i}, ...
                     'DimensionNames',{'time'},...
                     'Unit', unit{i}, ...
                     'LongVariableName',long_varname{i});
         end
end
