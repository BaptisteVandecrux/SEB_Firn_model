function []  = WriteVarToNC(namefile, data, varargin)
% This is an attempt to create a robust netcdf writting script
% It adds increment to the prescribed file name when file already exists or
% is not compatible with the variable that is being written (f.e. when
% dimensions do not match with existing file).
% The scripts can add Unit, LongVariableName and Comment as variable
% attribute.
%
% Baptiste Vandecrux bav@geus.dk
% ====================================================================


     flag = 0;
     count = 1;
     
     % set default values for the variable name, dimension size and names
     % and attributes
     param = struct;
     param.VariableName = 'var';
     param.Dimension = size(data);
     param.DimensionNames = {'dim1', 'dim2'};
     param.Unit = {};
     param.LongVariableName = {};
     param.Comment = {};
     
     if length(varargin)/2 ~= floor(length(varargin)/2)
         error('Wrong number of arguments')
     end
     
     for i = 1:2:length(varargin)
         if ~isfield(param,varargin{i})
             error('Wrong attribute name')
         else
             param.(varargin{i}) = varargin{i+1};
         end
     end
     for i = 1:length(param.DimensionNames)
         dim{2*i - 1} = param.DimensionNames{i};
         dim{2*i} = param.Dimension(i);
     end
         
         
     while flag == 0
        % nccreate can fail when a file already exist and has either 1) not
        % the appropriate definition for the variable dimensions or 2)
        % already a variable with that name in. In that situation we add an
        % increment to the file name.
        try
          nccreate(namefile,param.VariableName,...
          'Dimensions',dim,...
          'Format','classic') 
            flag = 1;
        catch
            count = count + 1;
            if count == 2
                namefile = strcat(namefile(1:(end-3)),sprintf('_%i.nc',count));
            elseif count < 10
                namefile = strcat(namefile(1:(end-5)),sprintf('_%i.nc',count));
            else
                namefile = strcat(namefile(1:(end-6)),sprintf('_%i.nc',count));
            end
        end
     end
 
    ncid = netcdf.open(namefile,'NC_WRITE'); %id of the nc file
    ID = netcdf.inqVarID(ncid,param.VariableName); %id of the variable in the file
    netcdf.putVar(ncid,ID,data);
    netcdf.close(ncid);

    ncid2 = netcdf.open(namefile,'NC_NOWRITE');
    ID = netcdf.inqVarID(ncid,param.VariableName);
    data_copy = netcdf.getVar(ncid2,ID);

    % prints a warning when one of the variabe contains nan
    if ~isequal(data,data_copy)
        fprintf('%s: check for NaNs\n',param.VariableName);
    end
    netcdf.close(ncid);
    if ~isempty(param.Unit)
        ncwriteatt(namefile,param.VariableName,'units',param.Unit);
    end
    if ~isempty(param.LongVariableName)
        ncwriteatt(namefile,param.VariableName,'long_name',param.LongVariableName);
    end
    if ~isempty(param.Comment)
        ncwriteatt(namefile,param.VariableName,'comment',param.Comment);
    end
end