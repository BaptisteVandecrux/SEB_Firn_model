function [tableout, subl] = ImportSnowpitData(c)


[~, ~, raw, dates] = xlsread( ['..\AWS_Processing\Input\Greenland_snow_pit_SWE.xlsx'],...
    'snow_pit_SWE_compiled_by_J_Box_','A2:P326','',@convertSpreadsheetExcelDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,12,13,14,16]);
raw = raw(:,[2,3,4,5,6,7,8,9,10,11]);
dates = dates(:,15);

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

data = reshape([raw{:}],size(raw));
tableout = table;
tableout.Station = cellVectors(:,1);
tableout.Date = datestr(datenum(data(:,6),data(:,5),data(:,4)));
tableout.SWE_pit = data(:,10);
clearvars data raw dates cellVectors R;

% Loading the sublimation estimates
filename = ['../AWS_Processing/Input/Sublimation estimates/', c.station, '_sublimation.txt'];

if exist(filename)==0
    disp('WARNING: Using sublimation estimate from NASA-E')
    filename = ['C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\'...
    'Code\AWS_Processing/Input/Sublimation estimates/NASA-E_sublimation.txt'];
end
delimiter = ';';
formatSpec = '%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
subl = table(dataArray{1:end-1}, 'VariableNames', {'time','estim'});
clearvars filename delimiter formatSpec fileID dataArray ans;
subl.time = datenum(subl.time,1,1); 

end
