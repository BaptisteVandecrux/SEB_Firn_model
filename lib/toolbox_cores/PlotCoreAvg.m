function [] = PlotCoreAvg(varargin) 
% PlotCore - Plot core data with different options
% PlotCore(SingleCore) 
%     where SingleCore is a a struct containing Info and Data from one specific core
% Plot(Core, KeyWord)
%     where Core is the core dataset and KeyWord is a string among:
%         'all' - Plots all cores in the dataset
%         'KAN-U' or 'KAN_U' - Plots all cores drilled at KAN-U site
%         'Summit' - Plots all cores drilled at Summit site
%         'EastGRIP' - Plots all cores drilled at EastGRIP site
%         'EKT' - Plots all cores drilled at EKT site
%         'Dye-2' - Plots all cores drilled at Dye-2 site
%         'Saddle' - Plots all cores drilled at Saddle site
%         'NASA-SE' - Plots all cores drilled at NASA-SE site
% Plot(Core, ind)
%     where is the Core dataset and ind an array of index for which core
%     will be plotted in windows of 5 subplots
% PlotCore(Core, KeyWord, ind) 
%     ##NOT RECOMMENDED##  Use "Plot(Core, ind)" to make subplots and use
%     CorrStudy to make overlap plots + differences in density.
%     where is the Core dataset and ind an array of index for which core
%     will be plotted, KeyWord can be either:
%         'subplot' to plot cores side by side (default setting 5 by
%         window)
%         'overlap' to overlap two core data given by ind 

    switch nargin
        case 3
            Core=varargin{1};
            KeyWord=varargin{2};
            cn=varargin{3};
        case 1
            Core{1}=varargin{1};
            KeyWord='subplot';
            cn=1;
        case 2
            if ~ischar(varargin{2})
                cn = varargin{2};
                Core = varargin{1};
                KeyWord='subplot';
            else
                KeyWord=varargin{2};
                Core=varargin{1};
            
                switch KeyWord
                    case 'all'
                        for i=0:floor(length(Core)/5)
                            cn=i*5 + (1:5);
                            PlotCoreAvg(Core,'subplot',cn);
                        end
                        return
                     
                    otherwise
                        if isnumeric(KeyWord)
                            year=KeyWord;
                            ind=[];
                            for i=1:length(Core)
                                if year == Core{i}.Info.DateCored.Year
                                    ind = [ind i];
                                end
                            end
                        else
                            name=KeyWord;
                            ind=[];
                            for i=1:length(Core)
                                info = Core{i}.Info.NearestCodeLocation;
                                if ~isempty(strfind(info,name))
                                    ind = [ind i];
                                end
                            end
                        end
                        if isempty(ind)
                                disp(sprintf(['KeyWord not recognized, please choose amongst: \n', ...
                                '''all'' - Plots all cores in the dataset \n',...
                                '''KAN-U'' or ''KAN_U'' - Plots all cores drilled at KAN-U site \n',...
                                '''Summit'' - Plots all cores drilled at Summit site\n',...
                                '''EastGRIP'' - Plots all cores drilled at EastGRIP site\n',...
                                '''EKT'' - Plots all cores drilled at EKT site\n',...
                                '''Dye-2'' - Plots all cores drilled at Dye-2 site\n',...
                                '''Saddle'' - Plots all cores drilled at Saddle site\n',...
                                '''NASA-SE'' - Plots all cores drilled at NASA-SE site\n',...
                                '''2015'' - Plots cores drilled in 2015 (same for other years)']))
                            end
                        temp = [];
                        for i = 1 : length(ind)
                            temp = [temp ind(i)];
                            if length(temp) == 5
                                PlotCoreAvg(Core,'subplot',temp);
                                temp = [];
                            end
                        end
                        if ~isempty(temp)
                            PlotCoreAvg(Core,'subplot',temp);
                        end
                        clear temp
                        return
                end
            end
        case 0
            PlotCore('all',Core);
            return
    otherwise
        disp('too many arguments for PlotCore');
    end
        
    aux=length(cn);
%     figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:(aux)
        if cn(i)<=length(Core)
            disp(sprintf('Plotting core %i',cn(i)))
            type=Core{cn(i)}.Data.Type;
            depth=Core{cn(i)}.Data.Depth;
            typeperc=Core{cn(i)}.Data.Type_perc;
            if isfield(Core{cn(i)}.Data,'Comment')
                comment=Core{cn(i)}.Data.Comment;
            end
            
%             subplot(1,aux,i)
            color = [0.7 0.8 1];
            edge1 = 0.01;
            edge2 = 1100;
            linecolor='b';
           
            hold on

% =================== Plotting stratigraphy ===========================
% when there is averaged stratigraphy
%             for j=1:length(depth)-1
%                 if strcmp(type{j},'ice') & ~isnan(depth(j)) & ~isnan(depth(j+1))                              
%                     if ~isnan(typeperc{j}) & isnumeric(typeperc{j})
%                         % if there is type_perc then scales the rectangle
%                         % to the percentage of ice in the layer
%                         rectangle('Position',...
%                             [edge1 depth(j) ...
%                             edge2*typeperc{j}/100 (depth(j+1)-depth(j))],...
%                             'FaceColor',color,'EdgeColor','none')
%                     else
%                         rectangle('Position',...
%                             [edge1 depth(j) ...
%                             edge2 (depth(j+1)-depth(j))],...
%                             'FaceColor',color,'EdgeColor','none')
%                     end
%                 end
%             end
            
            %Plotting bars showing the extent of the density profiles
            %used for average
            if Core{cn(i)}.Info.NumAvg > 1
                sorted = sort(Core{cn(i)}.Info.LengthAvg);
                for kk = 1:Core{cn(i)}.Info.NumAvg
                    rectangle('Position',...
                        [edge2-30-40 *(kk-1), 0,...
                        30 depth(sorted(kk))],...
                        'FaceColor',[0 0 0],'EdgeColor','none')
                end 
            end

            %Filling the rest of the plot window with grey
            aux2=depth(~isnan(depth));
            rectangle('Position',[0.01 aux2(end) edge2 (2000-aux2(end))],...
                'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')

            set(gca,'ydir','rev')
            title(sprintf('Name: %s \n Location: %s \n Nbr. of density profiles: %i.',...
                Core{cn(i)}.Info.Name, Core{cn(i)}.Info.NearestCodeLocation,...
                Core{cn(i)}.Info.NumAvg),'Interpreter', 'none')

            xlabel('Density (kg/m^3)')
            ylabel('Depth (cm)')
            ylim([0 2000])        
            xlim([0 edge2])
            
% ===================== Plotting Densities ===========================
            if strcmp(Core{cn(i)}.Info.Densities,'n')
                disp(sprintf('No density for Core %s',Core{cn(i)}.Info.Name));
            else
                density=MySmooth(Core{cn(i)}.Data.Density, 10);
                stdev=MySmooth(Core{cn(i)}.Data.DensStd,20);
                
                indend = min ([length(density),length(depth),length(stdev)]);
                plot(density(1:indend),depth(1:indend),'k','LineWidth',1.2)
                plot(density(1:indend)-stdev(1:indend),...
                    depth(1:indend),'r','LineWidth',0.5)
                plot(density(1:indend)+stdev(1:indend),...
                    depth(1:indend),'r','LineWidth',0.5)
            end
        else
            disp('Index exceeding dataset size')
        end               
    end
end