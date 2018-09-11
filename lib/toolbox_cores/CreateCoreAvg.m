function [CoreAvg] = CreateCoreAvg (Core)
% ======== Averaging cores of same year and same location =====
% Creating an new core register (CoreAvg) by averaging the duplicate cores
% from the old register (Core). The number of duplicate cores is stored in
% NumAvg (0 if no density profile, 1 if single core).

i= 1; % index in old register
count = 1; % index in new register
while i <= length(Core)
    ind = [];
    for j = 1: length(Core)
        % if it has density profile
        if ~strcmp(Core{i}.Info.Densities,'n')
            % if it has a specific place
            if ~isempty(Core{i}.Info.NearestCodeLocation)
                % if same year
                if Core{i}.Info.DateCored.Year == Core{j}.Info.DateCored.Year
                    % if same place
                    if strcmp(Core{i}.Info.NearestCodeLocation,...
                            Core{j}.Info.NearestCodeLocation  )
                        ind = [ind j];
                    end
                end
            end
        end
    end
    CoreAvg{count}= Core{i};
    
    dim = 0;
    count2 = 0;
    
    CoreAvg{count}.Info.LengthAvg = [];
    
    for j = ind
        if ~strcmp(Core{j}.Info.Densities,'n')
            %counts the number of duplicate cores that have density profiles
            count2 = count2 + 1;
            
            % max length of the duplicate density profiles
            temp = length(Core{j}.Data.Density);
            if temp > dim
                dim = temp;
                indmaxlength = j; %saving index of longest density profile
            end
            
            % saves the extent of each density profile used in the average
            for kk = length(Core{j}.Data.Depth):-1:1
                %goes backward in the density profile and saves the highest
                %index that is not nan
                    if ~isnan(Core{j}.Data.Density(kk))
                        indmax=kk;
                        break
                    end
            end
            CoreAvg{count}.Info.LengthAvg = ...
                [CoreAvg{count}.Info.LengthAvg indmax];
        end
    end
    clear temp

    if count2 > 1
        % if duplicate density profiles
        temp2 = NaN(max(CoreAvg{count}.Info.LengthAvg),count2);
        for j= 1:length(ind)
                temp2(1:length(Core{ind(j)}.Data.Density), j) = ...
                    Core{ind(j)}.Data.Density;
        end

        CoreAvg{count}.Data.Density = nanmean(temp2, 2);
        CoreAvg{count}.Data.DensStd = nanstd(temp2, 0, 2);
        clear temp2
        %replaces 0 standard deviation (i.e. sections where only one of the
        %replicate cores have density) by meas. uncertainty (50 kg/m^3)
        CoreAvg{count}.Data.DensStd(CoreAvg{count}.Data.DensStd==0)=50;
        CoreAvg{count}.Info.NumAvg = count2;
        if count==13
            disp(indmaxlength)
        end
        CoreAvg{count}.Data.Depth = Core{indmaxlength}.Data.Depth;
    else
      if strcmp(Core{i}.Info.Densities,'n')
        CoreAvg{count}.Data.DensStd = [];
        CoreAvg{count}.Info.Densities = 'n';
        CoreAvg{count}.Info.NumAvg = 0;
      else
        CoreAvg{count}.Data.DensStd = ...
            ones(length(Core{i}.Data.Density),1) * 50;
        CoreAvg{count}.Info.NumAvg = 1;
      end
    end

    CoreAvg{count}.Info.Name = sprintf('CoreAvg_%i_%i',...
        Core{i}.Info.DateCored.Year,count);
    clear temp2 
    
    fprintf('CoreAvg %i now contains',count)
    if length(ind) > 1
        disp(ind)
        i = ind(end) + 1;
    else
        disp(i)
        i = i+1;
    end
    count = count + 1;

end
end