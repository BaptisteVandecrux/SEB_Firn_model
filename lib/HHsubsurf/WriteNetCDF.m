function []  = WriteNetCDF(namefile, data, varname, dim1, nblayer, unit)
     flag = 0;
     count = 1;
     if nblayer == 1
         namevar2='var';
     else
         namevar2= 'layers';
     end
     
     while flag == 0
        try
          nccreate(namefile,varname,...
          'Dimensions',{namevar2,nblayer,'timesteps',dim1},...
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
    ID = netcdf.inqVarID(ncid,varname); %id of the variable in the file
    netcdf.putVar(ncid,ID,data);
    netcdf.close(ncid);

    ncid2 = netcdf.open(namefile,'NC_NOWRITE');
    ID = netcdf.inqVarID(ncid,varname);
    data_copy = netcdf.getVar(ncid2,ID);

    if ~isequal(data,data_copy)
        fprintf('%s: check for NaNs\n',varname);
    end
    netcdf.close(ncid);
    ncwriteatt(namefile,varname,'units',unit);
end