function [] = writetable2csv(data, filename)
        data_mat = table2array(data);
        data_mat(isnan(data_mat))=-999;

        for j = 1:size(data,2)
            data.(j) = data_mat(:,j);
        end
        varnames = data.Properties.VariableNames;
          fid = fopen(filename,'w');
          for i = 1:length(varnames)
            if i ==length(varnames)
                fprintf(fid, sprintf('%s\n',varnames{i}));
            else
                fprintf(fid, sprintf('%s; ',varnames{i}));
            end
          end

          fclose(fid);
          M = table2array(data);
          M(isnan(M))=-999;
        FastDLMwrite (filename, M, ';');
end
