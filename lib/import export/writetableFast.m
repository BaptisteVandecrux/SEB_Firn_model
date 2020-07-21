function [] = writetableFast(data, filename, delim)

    varnames = data.Properties.VariableNames;
	
	fid = fopen(filename,'w');
	for i = 1:length(varnames)
	    if i ==length(varnames)
	        fprintf(fid, sprintf('%s\n',varnames{i}));
	    else
	        fprintf(fid, sprintf('%s%s',varnames{i},delim));
	    end
	end
	
	fclose(fid);
	M=table2array(data);
	M(isnan(M)) = -999;
	FastDLMwrite(filename, M, delim);
end
    