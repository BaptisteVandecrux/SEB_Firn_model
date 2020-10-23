%      [STATUS,sheet_list,FORMAT] =xlsfinfo('C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Code\AWS_Processing\Input\ErrorData_AllStations.xlsx');
     [STATUS,sheet_list,FORMAT] =xlsfinfo('.\Input\ErrorData_AllStations.xlsx');
     
for i = 1:length(sheet_list)
%     [~, ~, raw] = xlsread('C:\Users\bav\OneDrive - Geological survey of Denmark and Greenland\Code\AWS_Processing\Input\ErrorData_AllStations.xlsx',...
    [~, ~, raw] = xlsread('.\Input\ErrorData_AllStations.xlsx',...
        sheet_list{i});
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,[2,3,1]);
    dates1 = cellVectors(:,1);
    dates2 = cellVectors(:,2);
    new_date1 = {};
    for ii = 1:length(dates1)
        try
            new_date1{ii} = datestr(datenum(dates1{ii,:},'dd-mm-yyyy HH:MM:SS'), 'yyyy-mm-dd HH:MM:SS');
        catch me
            try 
                new_date1{ii} = datestr(datenum(dates1{ii,:},'dd-mm-yyyy'), 'yyyy-mm-dd HH:MM:SS');
            catch me
                new_date1{ii} = datestr(datenum(dates1{ii,:},'dd-mmm-yyyy HH:MM:SS'), 'yyyy-mm-dd HH:MM:SS');      
            end
        end
    end
    new_date1=new_date1';
    new_date2 = {};
    for ii = 1:length(dates2)
        try
            new_date2{ii} = datestr(datenum(dates2{ii,:},'dd-mm-yyyy HH:MM:SS'),...
            'yyyy-mm-dd HH:MM:SS');
        catch me
            try 
                new_date2{ii} = datestr(datenum(dates2{ii,:},'dd-mm-yyyy'), 'yyyy-mm-dd HH:MM:SS');
            catch me
                new_date2{ii} = datestr(datenum(dates2{ii,:},'dd-mmm-yyyy HH:MM:SS'), 'yyyy-mm-dd HH:MM:SS');      
            end
        end
    end
    new_date2=new_date2';
    newCellVector = [new_date1, new_date2, {cellVectors{:,3}}', ...
        repmat({'UNKNOWN'},length(new_date2),1),...
        repmat({'Manually flagged by bav'},length(new_date2),1),...
        repmat({''},length(new_date2),1)];
    % Convert cell to a table and use first row as variable names
    T = cell2table(newCellVector,'VariableNames',...
        {'t0', 't1', 'variable', 'flag', 'comment', 'URL_graphic'});
    M = readtable([sheet_list{i}, '.csv']);
    M(:,end) = [];
    for k = 1:height(M)
        M.URL_graphic{k} = '';
    end
    % Write the table to a CSV file
    writetable([M;T],[sheet_list{i}, '.csv'])
end
