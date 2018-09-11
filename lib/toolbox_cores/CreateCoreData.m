function [Core] = CreateCoreData()
    % =================== Import overall data description ===================
    [~, ~, raw, dates] = xlsread('./Data/firn_cores_2009_2012_2013_2015_2016_temp.xlsx','firn_cores_2009_2012_2013_2015_','A2:S56','',@convertSpreadsheetExcelDates);
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,[1,2,6,7,10,14,15,16,17,18,19]);
    raw = raw(:,[3,4,5,8,9,13]);
    dates = dates(:,[11,12]);

    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dates); % Find non-numeric cells
    dates(R) = {NaN}; % Replace non-numeric Excel dates with NaN

    % Create output variable
    data = reshape([raw{:}],size(raw));

    % Allocate imported array to column variable names
    core = cellVectors(:,1);
    nearestcodelocation = cellVectors(:,2);
    N1 = data(:,1);
    W1 = data(:,2);
    Z1 = data(:,3);
    N2 = cellVectors(:,3);
    W2 = cellVectors(:,4);
    Z2 = data(:,4);
    depth1 = data(:,5);
    corer = cellVectors(:,5);
    datecored = datetime([dates{:,1}].', 'ConvertFrom', 'Excel');
    datelogged = datetime([dates{:,2}].', 'ConvertFrom', 'Excel');
    T2mCwhilelogging = data(:,6);
    densities = cellVectors(:,6);
    drilledby = cellVectors(:,7);
    loggedby1 = cellVectors(:,8);
    densitiesby2 = cellVectors(:,9);
    DO18 = cellVectors(:,10);
    remarks = cellVectors(:,11);

    % For code requiring serial dates (datenum) instead of datetime, uncomment
    % the following line(s) below to return the imported dates as datenum(s).

    % datecored=datenum(datecored);
    % datelogged=datenum(datelogged);
    % Clear temporary variables
    clearvars data raw dates cellVectors R;

    % =================== Creation of the dataset ===================

    for i=1:(length(core)-1)
        Core{i}.Info.Name=core{i};
        Core{i}.Info.NearestCodeLocation=nearestcodelocation{i};
        Core{i}.Info.Coord1=[N1(i) W1(i) Z1(i)];
        Core{i}.Info.Coord2={N2{i} W2{i} Z2(i)};
        Core{i}.Info.DepthMax=depth1(i);
        Core{i}.Info.Corer=corer{i};
        Core{i}.Info.DateCored=datecored(i);
        Core{i}.Info.DateLogged=datelogged(i);
        Core{i}.Info.T2mCwhilelogging=T2mCwhilelogging(i);
        Core{i}.Info.Densities=densities{i};
        Core{i}.Info.DrilledBy=drilledby{i};
        Core{i}.Info.LoggedBy=loggedby1{i};
        Core{i}.Info.DensityBy=densitiesby2{i};
        Core{i}.Info.DO18=DO18{i};
        Core{i}.Info.Remarks=remarks{i};


        % Import the data
        [num, txt, ~] = xlsread('./Data/firn_cores_2009_2012_2013_2015_2016_temp.xlsx',i+1);
        
        % Allocate imported array to column variable names
        ind = strcmp(txt(1,:),'depth');
        Core{i}.Data.Depth = num(:,ind);

        ind = strcmp(txt(1,:),'break');
        Core{i}.Data.Break = txt(2:end,ind);

        ind = strcmp(txt(1,:),'percol');
        Core{i}.Data.Percol = txt(2:end,ind);

        ind = strcmp(txt(1,:),'type');
        Core{i}.Data.Type = txt(2:end,ind);

        ind = strcmp(txt(1,:),'type_perc');
        Core{i}.Data.Type_perc = num(:,ind);
        Core{i}.Data.Type_perc(isnan(Core{i}.Data.Type_perc))=100;

        ind = strcmp(txt(1,:),'weight');
        Core{i}.Data.Weight = num(:,ind);

        ind = strcmp(txt(1,:),'weight_perc');
        Core{i}.Data.Weight_perc = num(:,ind);

        ind = strcmp(txt(1,:),'unrel_weight_perc');
        Core{i}.Data.Unrel_weight_perc = txt(2:end,ind);

        ind = strcmp(txt(1,:),'length');
        Core{i}.Data.Length1 = num(:,ind);

        ind = strcmp(txt(1,:),'length_std');
        Core{i}.Data.Length_std = num(:,ind);

        ind = strcmp(txt(1,:),'diameter');
        Core{i}.Data.Diameter = num(:,ind);

        ind = strcmp(txt(1,:),'diameter_std');
        Core{i}.Data.Diameter_std = num(:,ind);

        ind = strcmp(txt(1,:),'vol_std');
        Core{i}.Data.Vol_std = num(:,ind);

        ind = strcmp(txt(1,:),'density');
        Core{i}.Data.Density = num(:,ind);
        Core{i}.Data.Density(Core{i}.Data.Density==0) = NaN;
        Core{i}.Data.Density(find(strcmp(Core{i}.Data.Unrel_weight_perc,'?'))) = NaN;

        % Clear temporary variables
        clearvars num txt;
    end

    clearvars core nearestcodelocation N1 W1 Z1 N2 W2 Z2 depth1 corer  
    clearvars datelogged T2mCwhilelogging densities drilledby loggedby1 
    clearvars densitiesby2 DO18 remarks datecored i

end
