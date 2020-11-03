function maintenance = ImportMaintenanceFile(station)
    opts = spreadsheetImportOptions("NumVariables", 24);
    opts.Sheet = station;
    opts.DataRange = "A2:X30";
    opts.VariableNames = ["DateddmmyyyyHHMM", "reported", "SR1beforecm", "SR1aftercm", "SR2beforecm", "SR2aftercm", "T1beforecm", "T1aftercm", "T2beforecm", "T2aftercm", "W1beforecm", "W1aftercm", "W2beforecm", "W2aftercm", "NewDepth1m", "NewDepth2m", "NewDepth3m", "NewDepth4m", "NewDepth5m", "NewDepth6m", "NewDepth7m", "NewDepth8m", "NewDepth9m", "NewDepth10m"];
    opts.SelectedVariableNames = ["DateddmmyyyyHHMM", "reported", "SR1beforecm", "SR1aftercm", "SR2beforecm", "SR2aftercm", "T1beforecm", "T1aftercm", "T2beforecm", "T2aftercm", "W1beforecm", "W1aftercm", "W2beforecm", "W2aftercm", "NewDepth1m", "NewDepth2m", "NewDepth3m", "NewDepth4m", "NewDepth5m", "NewDepth6m", "NewDepth7m", "NewDepth8m", "NewDepth9m", "NewDepth10m"];
    opts.VariableTypes = ["string", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 2], "EmptyFieldRule", "auto");
    maintenance = readtable("C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Code\GEUS model\Input\maintenance.xlsx", opts, "UseExcel", false);
    clear opts
    maintenance.date = maintenance.T1aftercm*NaN;
    i=1;
    while ~isempty(maintenance.DateddmmyyyyHHMM{i})
        try maintenance.date(i) = datenum(maintenance.DateddmmyyyyHHMM{i}, 'dd-mm-yyyy HH:MM:SS');
        catch me
            maintenance.date(i) = datenum(maintenance.DateddmmyyyyHHMM{i}, 'dd-mm-yyyy');
        end
        i=i+1;
    end
    maintenance(i:end,:)=[];
end