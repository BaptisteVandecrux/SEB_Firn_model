function []  = WriteNetCDF2(namefile, data, varargin)
     flag = 0;
     count = 1;
     
     % set default values
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
         try dim{2*i - 1} = param.DimensionNames{i};
         catch me 
             sdjgh= 0;
         end
         dim{2*i} = param.Dimension(i);
     end
         
         
     while flag == 0
        try
          nccreate(namefile,param.VariableName,...
          'Dimensions',dim,...
          'Format','classic') 
            flag = 1;
                        info = ncinfo(namefile);

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