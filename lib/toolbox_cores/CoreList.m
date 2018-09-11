function [] = CoreList(Core)
%Print list of cores contained in dataset "Core" with index and names
disp('===== List of Cores ======')
disp('Ind     Name        Location      Date')
    for i=1:length(Core)
            fprintf('%i \t%s\t\t%s\t%s\n',...
                i, Core{i}.Info.Name, strtrim(Core{i}.Info.NearestCodeLocation),...
                datestr(datenum(Core{i}.Info.DateCored),'mmm-yyyy'));
    end
disp('===== End of List ======')
end