function index = SearchCore(Core,type,coresearch)

% function that looks for the index "index2 of a core "coresearch" in the
% core dataset "Core"
index=[];
    for i=1:length(Core)
        if isfield(Core{i}.Info,type)
            if iscell(Core{i}.Info.(type))
                Core{i}.Info.(type) = Core{i}.Info.(type){1};
            end
            if ~isempty(strfind(Core{i}.Info.(type),coresearch)) 
%                 if strcmp(Core{i}.Info.(type),coresearch)
%                 disp(sprintf('%s from %s is no. %i.',...
%                     Core{i}.Info.Name,Core{i}.Info.NearestCodeLocation,i));
                index=[index i];
            end
        end
    end
    if isempty(index)
        fprintf('The core "%s" was not found.\n',coresearch);
    end

%ordering the cores in chronological order
dates = zeros(size(index));
for i = 1:length(index)
    dates(i) = datenum(Core{index(i)}.Info.DateCored);
end
[~, i_ordered] = sort(dates);
index = index(i_ordered);

end