function [data_avg] = TempAvg(data)
    data_avg = table;
    data_avg.Year = table2array(data(:,1));
    data_avg.avg_year = mean(table2array(data(:,2:end)),2);
    data_avg.avg_JJA = mean(table2array(data(:,7:9)),2);
end
