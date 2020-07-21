function [] = WriteNC_2D(namefile, time, d2, data, d2_name,d2_unit, ...
    d2_long_name, varname, unit, long_varname)
% This is an attempt to create a robust netcdf writting script.
% It adds increment to the prescribed file name when file already exists or
% is not compatible with the variable that is being written (f.e. when
% dimensions do not match with existing file).
% The scripts can add Unit, LongVariableName and Comment as variable
% attribute.
% Example:
% time =datenum(2019,1,1,1:1000,0,0);
% scale_D2 = 0.1:0.1:100; % could be fixed measurements depth, height or station number
% data = {peaks(1000) exp(peaks(1000)) -peaks(1000)};
% varname = {'var1','var2','var3'};
% long_varname = {'long_var1','long_var2','long_var3'};
% unit = {'m' 's' '-'};
% 
% % note that the name, unit and description of the second dimension needs to
% % be provided
% WriteNC_2D('Example_2d.nc', time, scale_D2, data,...
%             'depth', 'm', 'Depth below the initial surface position',...
%             varname, unit, long_varname)
%
% Baptiste Vandecrux bav@geus.dk
% ====================================================================
if ~iscell(data)
    data = {data};
    varname = {varname};
    unit = {unit};
    long_varname = {long_varname};
end
    
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

     if or(size(d2,1) ==1,size(d2,2) ==1)
%          if there is one depth scale for all the columns then we use it
%          as dimension
        if isrow(d2)
            d2 = d2';
        end
         WriteVarToNC(namefile,d2, ...
             'VariableName',d2_name, ...
             'DimensionNames',{d2_name},...
             'Unit', d2_unit, ...
             'LongVariableName',d2_long_name);
         
         for i=1:length(data)
                 WriteVarToNC(namefile,data{i}, ...
                     'VariableName',varname{i}, ...
                     'DimensionNames',{d2_name,'time'},...
                     'Unit', unit{i}, ...
                     'LongVariableName',long_varname{i});
         end
     else
        % Otherwise we use "level" as dimension and add depth as a variable
        level = 1:size(d2,1);
        WriteVarToNC(namefile,level', ...
            'VariableName','level', ...
            'DimensionNames',{'level'},...
            'Unit', '', ...
            'LongVariableName','Model level');
         
        for i=1:length(data)
             WriteVarToNC(namefile,data{i}, ...
                 'VariableName',varname{i}, ...
                 'DimensionNames',{'level','time'},...
                 'Unit', unit{i}, ...
                 'LongVariableName',long_varname{i});
        end

        WriteVarToNC(namefile,d2, ...
            'VariableName',d2_name, ...
            'DimensionNames',{'level','time'},...
            'Unit', 'm', ...
            'LongVariableName','Depth below the surface');
     end

end