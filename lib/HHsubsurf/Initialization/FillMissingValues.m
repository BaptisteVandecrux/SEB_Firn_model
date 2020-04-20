function [missingnum, variable] = FillMissingValues(variable,method,opt)
    missing = find(isnan(variable(:,1)));
    if sum(missing) > 0
        notmissing = find(~isnan(variable(:,1)));
        if exist('opt')
        variable(missing) = interp1(notmissing,variable(notmissing),missing, method,opt);
        else
        variable(missing) = interp1(notmissing,variable(notmissing),missing, method);
        end
        missingnum = length(missing);
    else
        missingnum=0;
    end
end