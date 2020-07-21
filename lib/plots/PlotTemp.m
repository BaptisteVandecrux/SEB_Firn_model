function col = PlotTemp (varargin)
% function that allows to plot colored maps of any subsurface variable, few
% paired arguments are available to customize the result. More can be added
%
% Example:
%     figure
%     PlotTemp(time_mat, depth_mat_obs, mat_obs,...
%         'PlotTherm', 'yes',...
%         'PlotIsoTherm', 'yes',...
%         'IndNaN', [53324 35659]);
% or: 
%     f = figure('Visible', vis);
%     col = PlotTemp(TT,depth_act,rho_all,...
%     'PlotTherm', 'no',...
%     'PlotIsoTherm', 'yes',...
%     'ShowLegend','on',...
%     'cmap','jet',...
%     'XLabel','Time',...
%     'YLabel','Depth (m)',...
%     'CLabel','Density of frozen material (kg/m^3)',...
%     'Range', 350:25:900);
% 
%     col.TickLength=col.TickLength*5;
%     temp = col.YTickLabel;
%     for i = 1:length(temp)
%         if i/2==floor(i/2)
%             temp(i,:)='   ';
%         end
%     end
%     set(gca,'Color',[0.95 0.95 0.95]);

time_mat = varargin{1};
depth_mat = varargin{2};
mat = varargin{3};

if size(mat,1) ~= size(depth_mat,1)
    mat = vertcat(mat, mat(end,:));
    time_mat = vertcat(time_mat,time_mat(end,:));
end
    

% default parameters
param.PlotTherm = 'yes';
param.PlotIsoTherm = 'no';
param.XLabel = 'Time';
param.YLabel = 'Depth (m)';
param.Range = min(min(mat)):(max(max(mat))-min(min(mat)))/10:max(max(mat));
param.IndNaN = [];
param.ShowLegend = 'on';
param.cmap = 'jet';
param.CLabel = 'Observed Temperature (deg C)';
param.DepthSensor =  depth_mat;
param.DataIsoTherm = mat;

param.Tsurf = [];
param.ValueIsoTherm = [-1 -5];
param.Interp = 'off';
param.FlatSurface = 'no';

%assigning user defined parameters
if length(varargin)>3
    if floor((length(varargin)-3)/2) ~= (length(varargin)-3)/2
        fprintf('Error: Wrong number of arguments\n')
    end
    for i = 4:2:length(varargin)
        param.(varargin{i}) = varargin{i+1};
    end
end

if strcmp(param.FlatSurface,'yes')
    lm = fitlm(time_mat(1,:),depth_mat(1,:));
    trend_surf = lm.Coefficients.Estimate(2);
    fprintf('Average advection: %0.2f metres per decade\n',trend_surf*365*10)

    depth_mat =  depth_mat...
        -(trend_surf*time_mat+lm.Coefficients.Estimate(1));
end

% making sure that the data for isotherm have the right size
if size(param.DataIsoTherm,1) ~= size(depth_mat,1)
    param.DataIsoTherm = vertcat(param.DataIsoTherm, param.DataIsoTherm(end,:));
end

if ~isempty(param.Tsurf)
    depth_mat = [zeros(1,size(depth_mat,2)); depth_mat];
    time_mat = [time_mat(1,:); time_mat];
    mat = [param.Tsurf; mat];
    for i = 1 : size(depth_mat,2)
        indnan = or(isnan(depth_mat(:,i)), isnan(mat(:,i)));
        if sum(~indnan)>3
            depth_mat(indnan,i) = interp1(find(~indnan),depth_mat(~indnan,i),...
                find(indnan),'linear');
            mat(indnan,i) = interp1(find(~indnan), mat(~indnan,i),...
                find(indnan), 'linear');
        end
    end
end

if strcmp(param.Interp,'on')
    contourf(time_mat, depth_mat, mat, param.Range, 'Linestyle','none')
else
    surf(time_mat,depth_mat,zeros(size(mat)),mat,'EdgeColor','none');
%     pcolor(time_mat,depth_mat,mat)
    shading flat
end
grid off
set(gca,'layer','top')
view(0, -90)
axis tight
xlimit=get(gca,'XLim');
hold on
leg_text = cell(3,1);

if ~isempty(param.Tsurf)
    depth_mat = depth_mat(2:end,:);
    time_mat = time_mat(2:end,:);
    mat = mat(2:end,:);
end

if strcmp(param.PlotIsoTherm,'yes')
    if length(param.ValueIsoTherm)==1
        [~,h1]=contour(time_mat, depth_mat,param.DataIsoTherm,...
        [param.ValueIsoTherm(1) param.ValueIsoTherm(1)],'LineWidth',2,'Color','k');
        h2 = [];
        leg_text{1} = 'Pore close off density';
        leg_text{2} = '';
        
    elseif length(param.ValueIsoTherm)==2
        [~,h1]=contour(time_mat, depth_mat,param.DataIsoTherm,...
        [param.ValueIsoTherm(1) param.ValueIsoTherm(1)],'LineStyle',':','Color',[0.5 0.5 0.5]);
        [~,h2]=contour(time_mat, depth_mat,param.DataIsoTherm,...
        [param.ValueIsoTherm(2) param.ValueIsoTherm(2)],'LineStyle','-.','Color',[0.5 0.5 0.5]);
        leg_text{1} = '-1 ^o C isotherm';
        leg_text{2} = '-5 ^o C isotherm';
    end
end

if strcmp(param.PlotTherm,'yes')
    time_mat_2 = repmat(time_mat(1,:),size(param.DepthSensor,1)+1,1);
    sensor= repmat([0:size(param.DepthSensor,1)]',1,size(param.DepthSensor,2));
    depth_mat_2 = [zeros(1,size(time_mat,2)); param.DepthSensor];

%     depth_mat_2(isnan(depth_mat_2))=0;
%     depth_mat_2(depth_mat_2<=0)=NaN;
%     sensor(isnan(depth_mat_2))=NaN;
%     depth_mat_2(:,param.IndNaN)=NaN;
    [~,~]=contour(time_mat_2,depth_mat_2,sensor,1:size(mat,1),'LineStyle','-','Color',[0.5 0.5 0.5]);
    [~,h3]=contour(time_mat_2,depth_mat_2,sensor,[3 3],'LineStyle','-','Color',[0.5 0.5 0.5]);
    leg_text{3} = 'Thermistor';
end

if strcmp(param.ShowLegend,'yes')
    if ~strcmp(param.PlotTherm,'yes') && strcmp(param.PlotIsoTherm,'yes')

        leg = legend([h1, h2] ,leg_text{1},leg_text{2}, 'Location','SouthWest');
        
    elseif strcmp(param.PlotTherm,'yes') && strcmp(param.PlotIsoTherm,'yes')
        leg = legend([h1 h2 h3],leg_text{1},leg_text{2},leg_text{3}, 'Location','SouthWest');
    elseif strcmp(param.PlotTherm,'yes') && ~strcmp(param.PlotIsoTherm,'yes')
        leg = legend([h3],leg_text{3}, 'Location','SouthWest');
    end
    leg.FontSize = 7.8;
    p = leg.Position;
    leg.Position = [p(1) p(2) .85*p(3) .6*p(4)];
end

col =contourcmap(param.cmap,param.Range,'colorbar','on');

view(0,-90)

time_mod = datenum(0,0,time_mat(1,:));

set_monthly_tick(time_mod);
box on
axis tight

xlabel(param.XLabel)
ylabel(param.YLabel)
ylabel(col,param.CLabel,'Interpreter','tex')

set (gca, 'YLim', [0, 12], 'XLim',xlimit,...
    'XMinorTick','on','YMinorTick','on');

end