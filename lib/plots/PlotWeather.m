function [] = PlotWeather (data, varargin)
% function plotting the following variables contained in the table 'data':
%   'T_2m_degC'
%   'RH_2m_perc',...
%   'AirTemperature2mC'
%   'RelativeHumidity2mPerc'
%   'AirPressurehPa'
%   'WindSpeed10mms'
% 
%   accept pair arguments such as:
%   PlotWeather (data, 'XLim', [x_1 x_2])
%       will limit the x axis of the subplot to [x_1 x_2]. Default: all x
%       axis is shown.
%   PlotWeather (data, 'vis','off')
%       will not display the plot windows. They will still be printed as
%       files. Default: 'on'
%   PlotWeather (data, 'OutputFolder',FolderName)
%       will print the plots in an existing directory specified as string in
%       FolderName (f.e. FolderName = './folder'; ). Default: plots are
%       saved in the current working directory.

param.NameFile = 'WeatherData';
param.XLim = [data.time(1) data.time(end)];
param.vis = 'on';
param.OutputFolder ='.';
param.VarList = {'ShortwaveRadiationDownWm2','ShortwaveRadiationUpWm2',...
    'AirTemperature2C','RelativeHumidity2Perc','AirPressurehPa',...
    'WindSpeed2ms','LongwaveRadiationDownWm2','LongwaveRadiationUpWm2',...
    'Snowfallmweq','Rainfallmweq'};
param.LabelList = {'Downward SW \newline radiation (W/m^2)',...
    '  Upward SW \newline radiation (W/m^2)',...
    'Air temperature \newline          (^oC)',...
    '    Relative \newline humidity (%)',...
    'Air Pressure \newline      (hPa)',...
    'Wind speed \newline      (m/s)',...
    'Downward LW \newline radiation (W/m^2)',...
    'Upward LW \newline radiation (W/m^2)',...
     'Snowfall rate \newline (m weq/hr)',...
     'Rainfall rate \newline (m weq/hr)'};
param.Type ='line';

for i = 1:2:length (varargin)
    param.(varargin{i}) = varargin{i+1};
end

if isempty(strfind([varargin{1:2:end}],'Origins'))
     for i = 1:length(param.VarList)
        param.Origins{i} = [param.VarList{i} '_Origin'];
     end
end

if strcmp(param.Type,'line')
    plot_f = @(x,y,varargin) plot(x,y,varargin{:});
elseif strcmp(param.Type,'stairs')
    plot_f = @(x,y,varargin) stairs(x,y,varargin{:});
end

i_remove = [];
for i = 1:length(param.VarList)
    if ~ismember(param.VarList{i},data.Properties.VariableNames)
        % not plotting variables in VarList that are not in data
        i_remove = [i_remove i];
    else
        % not plotting variables that are filled with nan
        if sum(~isnan(data.(param.VarList{i})))<10
            i_remove = [i_remove i];
        end
    end
end
param.VarList(i_remove) = [];
param.Origins(i_remove) = [];
param.LabelList(i_remove) = [];

%% Setting variable names and labels
all_levels = [];
for i = 1:length(param.Origins)
    if ismember(param.Origins{i},data.Properties.VariableNames)
        all_levels =[all_levels; unique(data.(param.Origins{i}))];
        if sum(all_levels<0) >0
            error(sprintf('Non valid origin code for:\n %s',param.Origins{i}))
        end
    else
            data.(param.Origins{i}) = zeros(size(data.time));
            all_levels =[all_levels; unique(data.(param.Origins{i}))];
    end
end
all_levels = sort(unique(all_levels));
all_levels = all_levels(~isnan(all_levels));
  
legend_text = {'Original data','adjusted from CP2',...
                    'adjusted from Swiss Camp','adjusted from KAN-U',...
                    'adjusted from RCM','calculated using MODIS albedo',...
                    'taken from previous year',...
                    'Reconstructed by Charalampidis et al. 2015',...
                    'Kobblefjord station','adjusted from NOAA station',...
                    'adjusted from Miller et al. 2017'};


%% Plotting

f = figure('Visible',param.vis);
num_plot = length(param.VarList);
ha = tight_subplot(num_plot,1,0.01, .07, 0.08);

count = 1;
for i = 1:length(param.VarList)
        if count <= num_plot
            set(f,'CurrentAxes',ha(count))
        else

            legend('additional data','shifted in time','old data')
            i_file = 1;
            NameFile = sprintf('%s/%s_%i.pdf',param.OutputFolder, param.NameFile, i_file);
            while exist(NameFile, 'file') == 2
                i_file = i_file + 1;
                NameFile = sprintf('%s/%s_%i.pdf',param.OutputFolder, param.NameFile,i_file)  ;
            end
            print(f,NameFile,'-dpdf');
            f = figure('Visible',param.vis);
            ha = tight_subplot(num_plot,1,0.01, [.07 .03], 0.05);
            count = 1;
            set(f,'CurrentAxes',ha(count))
        end
        
        set(gca,'fontsize',13)
        hold on
        

        for ii =1:length(all_levels)
            data2{ii} = data;
            data2{ii}.(param.VarList{i})(data.(param.Origins{i})~=all_levels(ii)) = NaN;
        end
        
        plot_f(data.time,data.(param.VarList{i}),'k') 
        for j = 1:length(all_levels)
            switch all_levels(j)
%                 case 0
%                 plot_f(data.time,data2{j}.(param.VarList{i}),...
%                     'k')                
                case 1
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'b')
                case 2
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'r')
                case 3
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'Color',[0,100,0]/255)          %RGB('dark green'))
                case 4
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'Color',[0.4940    0.1840    0.5560])
                case 5
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'Color',[255,165,0]/255)        %RGB('orange'))
                case 6
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'Color',[0,1,1])                %RGB('cyan'))
                case 7
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'Color',[139,69,19]/255)         %RGB('brown'))
                case 8
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'm')         %RGB('brown'))
                case 9
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'c')
                case 10
                plot_f(data.time,data2{j}.(param.VarList{i}),...
                    'g')
%             otherwise
%                 plot_f(data.time,data2{j}.(param.VarList{i}))
            end
        end

        axis tight

        if count/2 ==floor(count/2)
            set(gca,'YAxisLocation','right')
        end
        if count == 1
            % in the first subplot we set up the legend for the whole
            % figure
           
            ff = [];
            ylimits = get(gca,'YLim');
            xlimits = get(gca,'XLim');
            
            for jj = 1:length(all_levels)
                switch all_levels(jj)
                case 0
                ff = [ff, plot(1:2,1:2,'k')];                
                case 1
                ff = [ff, plot(1:2,1:2,'b')];
                case 2
                ff = [ff, plot(1:2,1:2,'r')];
                case 3
                ff = [ff, plot(1:2,1:2,'Color',[0,100,0]/255)];
                case 4
                ff = [ff, plot(1:2,1:2,'Color',[0.4940    0.1840    0.5560])];
                case 5
                ff = [ff, plot(1:2,1:2,'Color',[255,165,0]/255)];
                case 6
                ff = [ff, plot(1:2,1:2,'Color','c')];               
                case 7
                ff = [ff, plot(1:2,1:2,'Color',[139,69,19]/255)];
                case 8
                ff = [ff, plot(1:2,1:2,'m')];
                case 9
                ff = [ff, plot(1:2,1:2,'c')];
                case 10
                ff = [ff, plot(1:2,1:2,'g')];
                end
            end
            
            xlim(xlimits);
            ylim(ylimits);
            if length(ff)>1
                legendflex(ff, legend_text(all_levels+1), 'ref', gcf, ...
                       'anchor',  [2 6] , ...
                       'buffer',[0 -40], ...
                       'ncol',3, ...
                       'fontsize',13);
            else
                legendflex(legend_text(all_levels+1), 'ref', gcf, ...
                       'anchor',  [2 6] , ...
                       'buffer',[0 -40], ...
                       'ncol',3, ...
                       'fontsize',13);
            end
        end
        ylabel(param.LabelList{i},...
            'interpreter','tex',...
            'HorizontalAlignment','center')
%         if length(data.time)>365*24
            set_monthly_tick(data.time);
%         end
        if count<num_plot
            set(gca,'XTickLabel',[])
        end
        box on
        xlim(param.XLim);
        count = count +1;
end
for i = count:num_plot
    set(f,'CurrentAxes',ha(i))
    set(gca,'Visible','off')
end

i_file = 1;
    NameFile = sprintf('%s/%s_%i.pdf',param.OutputFolder, param.NameFile, i_file);
    while exist(NameFile, 'file') == 2
        i_file = i_file + 1;
        NameFile = sprintf('%s/%s_%i.pdf',param.OutputFolder, param.NameFile,i_file)  ;
    end
orient(f,'landscape')
print(f,NameFile,'-dpdf');    
end