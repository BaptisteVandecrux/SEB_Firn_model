function [f] = PlotCore(Core, varargin)
% PlotCore - Plot core data with different options
% PlotCore(SingleCore)
%     where SingleCore is a a struct containing Info and Data from one specific core
% Plot(Core, param.KeyWord)
%     where Core is the core dataset and param.KeyWord is a string among:
%         'all' - Plots all cores in the dataset
%         'KAN-U' - Plots all cores drilled at KAN-U site
%         'Summit' - Plots all cores drilled at Summit site
%         'EastGRIP' - Plots all cores drilled at EastGRIP site
%         'EKT' - Plots all cores drilled at EKT site
%         'Dye-2' - Plots all cores drilled at Dye-2 site
%         'Saddle' - Plots all cores drilled at Saddle site
%         'NASA-SE' - Plots all cores drilled at NASA-SE site
% Plot(Core, ind)
%     where is the Core dataset and ind an array of index for which core
%     will be plotted in windows of 5 subplots
% PlotCore(Core, param.KeyWord, ind)
%     ##NOT RECOMMENDED##  Use "Plot(Core, ind)" to make subplots and use
%     CorrStudy to make overlap plots + differences in density.
%     where is the Core dataset and ind an array of index for which core
%     will be plotted, param.KeyWord can be either:
%         'subplot' to plot cores side by side (default setting 5 by
%         window)
%         'overlap' to overlap two core data given by ind

%% default settings:
param.PlotStrat = 'yes';
param.Ylim = 20;
param.vis = 'on';
param.CoreNumber = 1:length(Core);
param.Site = [];
param.PlotLineH = [];
param.Year = [];
param.Order = '';
param.Color ='k';
% Assigning the user's option to the option structure (param)
for i = 1:2:length (varargin)
    param.(varargin{i}) = varargin{i+1};
end

%% Selection by year
% If a year is is specified then indexes of the cores from that years are
% found.
if ~isempty(param.Year)
    param.CoreNumber = [];
    for i=1:length(Core)
        if param.Year == Core{i}.Info.DateCored.Year
            param.CoreNumber = [param.CoreNumber i];
        end
    end
end

%% Selection by site
if ~isempty(param.Site)
    if strcmp(param.Site, 'all')
        for i=1:floor(length(Core)/5)+1
            f = PlotCore(Core,'CoreNumber',(i-1)*5 + (1:5));
        end
        return
    else
        param.CoreNumber = [];
        for i=1:length(Core)
            info = Core{i}.Info.NearestCodeLocation;
            if ~isempty(strfind(info,param.Site))
                param.CoreNumber = [param.CoreNumber i];
            end
        end
    end
end

%% Selection by core number
if isempty(param.CoreNumber)
    disp(sprintf(['param.KeyWord not recognized, please choose amongst: \n', ...
        '''all'' - Plots all cores in the dataset \n',...
        '''KAN-U'' - Plots all cores drilled at KAN-U site \n',...
        '''Summit'' - Plots all cores drilled at Summit site\n',...
        '''EastGRIP'' - Plots all cores drilled at EastGRIP site\n',...
        '''EKT'' - Plots all cores drilled at EKT site\n',...
        '''Dye-2'' - Plots all cores drilled at Dye-2 site\n',...
        '''Saddle'' - Plots all cores drilled at Saddle site\n',...
        '''NASA-SE'' - Plots all cores drilled at NASA-SE site\n',...
        '''Crawford Point'' - Plots all cores drilled at NASA-SE site\n',...
        '2015 - Plots cores drilled in 2015 (same for other years)']))
    return
end

%% ordering the cores in prescribed order
switch param.Order
    case 'Chronological'
    dates = zeros(size(param.CoreNumber));
    for i = 1:length(param.CoreNumber)
        dates(i) = datenum(Core{param.CoreNumber(i)}.Info.DateCored);
    end
    [~, i_ordered] = sort(dates);
    param.CoreNumber = param.CoreNumber(i_ordered);
    
    case 'Elevation'
    elev = zeros(size(param.CoreNumber));
    for i = 1:length(param.CoreNumber)
        elev(i) = datenum(Core{param.CoreNumber(i)}.Info.Coord1(3));
    end
    [~, i_ordered] = sort(elev);
    param.CoreNumber = param.CoreNumber(i_ordered);
end

%% Storing window in array
% if the function creates multiple windows then each one of them is
% stored in a cell array f
% temp = [];
% for i = 1 : length(param.CoreNumber)
%     temp = [temp param.CoreNumber(i)];
%     if length(temp) == 5
%         if exist('f','var')
%             f{length(f)+1} =  PlotCore(Core,'CoreNumber',temp);
%         else
%             f{1} = PlotCore(Core,'CoreNumber',temp);
%         end
%         temp = [];
%     end
% end
% temp(length(temp)+1:5) = length(Core)+1;
% if exist('f','var')
%     f{length(f)+1} =  PlotCore(Core,'CoreNumber',temp);
% else
%     f{1} = PlotCore(Core,'CoreNumber',temp);
% end
% clear temp
% return

%% Plotting
NumCores=length(param.CoreNumber);
if NumCores>8
    f = figure('Visible',param.vis,'units','normalized','outerposition',[0 0 1 1]) ;
else
    f = figure('Visible',param.vis,'units','normalized','outerposition',[0 0 NumCores/8+0.05 1]) ;
end
if NumCores>=15
    [ha,~] = tight_subplot(1,NumCores,0,[.13 .2],[0.05 .08]);
elseif NumCores<=2
    [ha,~] = tight_subplot(1,NumCores,0.003,[.13 .08],0.15);
else
    [ha,~] = tight_subplot(1,NumCores,0.003,[.13 .08],0.08);
end
hold on
for i=1:(NumCores)
    if i<=length(param.CoreNumber)
    if param.CoreNumber(i)<=length(Core)
        fprintf('Plotting core %i\n',param.CoreNumber(i));
        type=Core{param.CoreNumber(i)}.Data.Type;
        depth=Core{param.CoreNumber(i)}.Data.Depth/100;
        typeperc=Core{param.CoreNumber(i)}.Data.Type_perc;
        if isfield(Core{param.CoreNumber(i)}.Data,'Comment')
            comment=Core{param.CoreNumber(i)}.Data.Comment;
        end
        
        set(f,'CurrentAxes',ha(i))
        color = [0.7 0.8 1];
        edge1 = 0.01;
        edge2 = 950;
        linecolor='b';
        box on
        hold on
        
        %% =================== Plotting stratigraphy =======================
        if strcmp(param.PlotStrat,'on') || strcmp(param.PlotStrat,'yes')
            if ~isempty(type)
                for j=1:min(length(depth),length(type))-1
                    
                    if ~isnan(depth(j)) && ~isnan(depth(j+1))
                        frac = 0;
                        if ~isempty(type{j})
                            switch type{j}
                                case {'dry snow','snow'}
                                    color = RGB('cream');
                                case 'wet snow'
                                    color = RGB('light green');
                                case {'depth hoar'}
                                    color = RGB('light light green');
                                case 'unwetted firn'
                                    color = RGB('light blue');
                                case 'firn'
                                    color = RGB('blue');
                                case 'wetted firn'
                                    color = RGB('dark blue');
                                case 'ice'
                                    color = 'w';
                                case 'na'
                                    color = 'w';
                                otherwise
                                    color = 'w';
                                    fprintf('%s not recognized in color code\n',type{j});
                            end
                        else
                            color = 'w';
                        end
                        rectangle('Position',[edge1 depth(j) edge2 (depth(j+1)-depth(j))],...
                            'FaceColor',color,'EdgeColor','none');

                        % if ice, overwrites previous rectangle
                        if ~isnan(typeperc(j)) && isnumeric(typeperc(j))
                            frac = typeperc(j)/100;
                        end
                        if frac > 0
                            rectangle('Position',[edge1 depth(j) frac*edge2 (depth(j+1)-depth(j))],...
                                'FaceColor','k','EdgeColor','none');
                        end
                        
%                         if isfield(Core{param.CoreNumber(i)}.Data,'Comment')
%                             if ~isempty(strfind(comment{j},'ipe'))&& ~isnan(depth(j)) && ~isnan(depth(j+1))
%                                 rectangle('Position',...
%                                     [edge2-20 depth(j) edge2-10 (depth(j+1)-depth(j))],...
%                                     'FaceColor','r','EdgeColor','none')
%                             end
%                         end
                    end
                end
            end
        elseif strcmp(param.PlotStrat,'simple')
            if ~isempty(type)
                for j=1:min(length(depth),length(typeperc))-1
                    
                    if ~isnan(depth(j)) && ~isnan(depth(j+1))
                        if ~isnan(typeperc(j)) && isnumeric(typeperc(j))
                            frac = typeperc(j)/100;
                        else
                            frac = 0;
                        end
                        if frac > 0
                            rectangle('Position',[edge1 depth(j) frac*edge2 (depth(j+1)-depth(j))],...
                                'FaceColor',RGB('light blue'),'EdgeColor','none');
                        end
                    end
                end
            end
        end
        
        %% ================ Plotting density ===========================
        if strcmp(Core{param.CoreNumber(i)}.Info.Densities,'n')
            fprintf('No density for Core %s\n',...
                Core{param.CoreNumber(i)}.Info.Name);
        else
            density=Core{param.CoreNumber(i)}.Data.Density;
            density(density==0)=NaN;
            plot(density,depth(1:length(density)),param.Color,'LineWidth',1.7)

            plot([917 917],[0 max(depth(1:length(density)))],'--')
            if isfield(Core{param.CoreNumber(i)}.Data,'Density_origin')
                if sum(Core{param.CoreNumber(i)}.Data.Density_origin)>2
                    density_2 = density;
                    density_2(~Core{param.CoreNumber(i)}.Data.Density_origin) = NaN;
                    plot(density_2, ...
                        depth(1:length(density_2)), ...
                        'm','LineWidth',1.7)
                end
            end

        end
        
        %% ============= Plotting horizontal lines =====================
        if ~isempty(param.PlotLineH)
            if strcmp(Core{param.CoreNumber(i)}.Info.Densities,'y')
                % code for when things were in m weq
%                 depth_m_weq = cumsum(density*0.01)/1000;
%                 ind = find(~isnan(depth_m_weq),1,'last');
%                 depth_m_weq = depth_m_weq(1:ind);
% 
                numLine = length(param.PlotLineH)-1;
                col = hsv(max(3,numLine));

                if max(depth) >= param.PlotLineH(1)
                    for ii = 2:numLine+1
                        if max(depth) >= param.PlotLineH(ii)
                            [~, ind] = min(abs(depth-param.PlotLineH(1)));
                            hh(ii-1) = plot([0 900],depth([ind ind])+(ii-1)/10, ...
                                'LineWidth',3, 'Color',col(ii-1,:));
                        
                            [~, ind2] = min(abs(depth-param.PlotLineH(ii)));
                            plot([0 900],depth([ind2 ind2]), ...
                                'LineWidth',3, 'Color',col(ii-1,:))
                        else
                            sprintf('Core %s depth max %0.2f m \n lower bound of comparison is %0.2f\n',...
                                Core{param.CoreNumber(i)}.Info.Name,max(depth),param.PlotLineH(ii))
                        end
                    end
                else
                    sprintf('Core %s depth max %0.2f m \n upper bound of comparison is %0.2f\n',...
                            Core{param.CoreNumber(i)}.Info.Name, ...
                            max(depth), ...
                            param.PlotLineH(1))
                end
            end

        end            
        
        %% Plotting grey box from end of data until 200m deep
        if isnan(depth(end))
            aux2=min(depth(~isnan(depth)),depth(end));
        else
            aux2=depth(end);
        end
%         rectangle('Position',[0.02 aux2(end) 940 19.90],...
%             'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
        rectangle('Position',[0.01 aux2(end) 950 0],...
            'FaceColor','none','EdgeColor','k');

        %% Axis and title
        set(gca,'ydir','rev')
        
        switch param.Order
            case 'Chronological'
                title_text = sprintf('%s (%i)',Core{param.CoreNumber(i)}.Info.Name,...
            Core{param.CoreNumber(i)}.Info.DateCored.Year);
            case 'Elevation'
                title_text = sprintf('%s (%im)',Core{param.CoreNumber(i)}.Info.Name,...
            floor(Core{param.CoreNumber(i)}.Info.Coord1(3)));
            otherwise
                title_text = Core{param.CoreNumber(i)}.Info.Name;
        end
        h_title = title(title_text, 'Interpreter', 'none');
        h_title.Units = 'normalized';
        if NumCores>=15
            h_title.HorizontalAlignment ='left';
             set(h_title,'Rotation',45);
        elseif NumCores>=10
            if i/2==floor(i/2)
                h_title.Position = h_title.Position + [0 0.03 0];
            end
        else
            if i/2==floor(i/2)
                h_title.Position = h_title.Position + [0 0.04 0];
            end
        end

%         title(sprintf('Name: %s \n Location: %s \n Elevation: %i m',...
%             Core{param.CoreNumber(i)}.Info.Name,...
%             Core{param.CoreNumber(i)}.Info.NearestCodeLocation,...
%             Core{param.CoreNumber(i)}.Info.Coord1(3)),...
%             'Interpreter', 'none');
        if i==1
            xlabel('Density (kg/m^3)',...
                'Interpreter','tex',...
                'Position',[500*length(param.CoreNumber) 20.75]);
            ylabel('Depth (m)');
        
        elseif i == NumCores
            set(gca,'YAxisLocation','right')
        else
            set(gca,'YTickLabel','');
        end
        xlim([0 950]);
        if ~isempty(param.Ylim)
            set(gca,'YLim',[0 param.Ylim]);
        else
            set(gca,'YLim',[0 max(depth)]);
        end
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',12,...
            'layer','top')
    if NumCores>=10
        h_title.FontSize = 9;
    end
    else
        set(f,'CurrentAxes',ha(i))
        set(gca,'Visible','off')
        disp('Index exceeding dataset size')
    end
    else
        set(f,'CurrentAxes',ha(i))
        set(gca,'Visible','off')
    end
end
end