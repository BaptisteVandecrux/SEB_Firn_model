function [r, depth, densdiff] = OverlapPlot(Core, cn, varargin)
param.lag = 0;
param.span = 5;
param.PlotDifference = 'yes';
param.YLimit = [];
param.PlotStrat = 'yes';

% Assigning the user's option to the option structure (param)
for i = 1:2:length (varargin)
    param.(varargin{i}) = varargin{i+1};
end
set(0,'DefaulttextInterpreter','none')
if length(cn)~=2
    disp('Enter only 2 indexes.')
    return
end

if strcmp(param.PlotDifference,'yes')    
    ha = tight_subplot(1,2,0.02,0.1,0.1);
    axes(ha(1));
end
hold on

for i=1:2
        if i == 1
            type=Core{cn(i)}.Data.Type;
            typeperc=Core{cn(i)}.Data.Type_perc;
            depth=Core{cn(i)}.Data.Depth/100;
        else
            type=Core{cn(i)}.Data.Type(1:(end));
            typeperc=Core{cn(i)}.Data.Type_perc(1:(end));
            if param.lag>=0
                depth=Core{cn(i)}.Data.Depth((1+param.lag):end)/100;
            else
                depth=[(param.lag:0)'; Core{cn(i)}.Data.Depth]/100;
            end
        end           

        if i==1
            color = RGB('light light green');
            edge1=1;
        else
            color = [0.85 0.8 0.95];
            edge1=475;
        end
            width_box=475;

        %% Plotting stratigraphy
        if ~strcmp(param.PlotStrat,'no') 
        for j=1:min(length(depth),length(type))-1
            frac = 0;
            % if ice, overwrites previous rectangle
            if ~isnan(typeperc(j)) && isnumeric(typeperc(j))
                frac = typeperc(j)/100;
            end
            if frac > 0
                rectangle('Position',[edge1 depth(j) frac*width_box (depth(j+1)-depth(j))],...
                    'FaceColor',color,'EdgeColor','none');
            end
%             if i==1
%                 if depth(j)==1
%                     rectangle('Position',[edge1 depth(j)-3 width_box 3],...
%                         'FaceColor', 'k','EdgeColor','none')
%                 end
%             else
%                 if depth(j)==1
%                     rectangle('Position',[edge1 param.lag width_box 3],...
%                     'FaceColor', 'k','EdgeColor','none')
%                 end
%             end
        end
        end
end
            
        %% Plotting density
%         depth1=Core{cn(1)}.Data.Depth;
%         depth2 = Core{cn(2)}.Data.Depth;
%         thick1 = [depth2(1); depth1(2:end)-depth1(1:end-1)];
%         thick2 = [depth2(1); depth2(2:end)-depth2(1:end-1)];
%  i=1;
%  j=1;
%  common_depth_scale = [];
% while and(i<=length(thick1),j<=length(thick2))
%     if thick1(i)>=thick2(j)
%         common_depth_scale = [common_depth_scale, depth1(i)];
%         j = find(depth2>depth1(i));
%         if j>1
%             if depth2(j-1)<depth1(i)
%                 j = j+1;
%             end
%         end
%         i = i+1;
%     else
%         common_depth_scale = [common_depth_scale, depth2(j)];
%         i = find(depth1>depth2(j));
%         if i>1
%             if depth1(i-1)<depth2(j)
%                 i = i+1;
%             end
%         end
%         j = j+1;
%     end
% end
    
for i=1:2
        if i == 1
            depth=Core{cn(i)}.Data.Depth/100;
        else
            if param.lag>=0
                depth=Core{cn(i)}.Data.Depth((1+param.lag):end)/100;
            else
                depth=[(param.lag:0)'; Core{cn(i)}.Data.Depth]/100;
            end
        end     
        if i==1
            edge1=1;
            width_box=475;
            linecolor = RGB('dark green');
            posbox=[0.2 0.87 0.3 0.08];
        else
            edge1=475;
            width_box=950;
            linecolor=[0.3 0.2 0.35];
            posbox=[0.5 0.87 0.3 0.08];
        end
        
        if strcmp(Core{cn(i)}.Info.Densities,'n')
            fprintf('No density for Core %s\n',Core{cn(i)}.Info.Name);
        else
            density=Core{cn(i)}.Data.Density;
            density= MySmooth(density,param.span);
            temp=min(length(depth),length(density));
            hh(i) = plot(density(1:temp),depth(1:temp),'Color',linecolor,'LineWidth',1.7);
        end
end
%         grey box for no data
% [~, ind_min] = min([max(Core{cn(1)}.Data.Depth),max(Core{cn(2)}.Data.Depth)]);
% depth=Core{cn(ind_min)}.Data.Depth/100;
% aux2=depth(~isnan(depth));
% rectangle('Position',[1 aux2(end) width_box*2 max((20-aux2(end)),0)],...
%     'FaceColor',[0.85 0.85 0.85],'EdgeColor','none')
        
set(gca,'ydir','rev','layer','top')
if isempty(param.YLimit)
    ylim([min(param.lag,0) min(max(Core{cn(1)}.Data.Depth), max(Core{cn(2)}.Data.Depth))]) 
else
    ylim([min(param.lag,0), param.YLimit])
end
xlim([0 950])
for i =1:2
    TF = isstrprop(Core{cn(i)}.Info.Name,'digit');
    if sum(TF)<4
        Core{cn(i)}.Info.Name = strcat(Core{cn(i)}.Info.Name,'_',num2str(Core{cn(i)}.Info.DateCored.Year));
    end
end
hleg = legend(hh, ...
    {Core{cn(1)}.Info.Name, Core{cn(2)}.Info.Name},...
    'Location','NorthOutside');
legend('boxoff');
set(hleg,'interpreter','none')
set(gca,'XMinorTick','on','YMinorTick','on')
box on
ylimit = get(gca,'YLim');
ylim([0 min(length(Core{cn(1)}.Data.Density),...
    length(Core{cn(2)}.Data.Density))]/100)

    %% Plotting density difference
    if strcmp(Core{cn(i)}.Info.Densities,'n')
        depth=NaN;
        densdiff=NaN;
        r=NaN;
    else
        d1=Core{cn(1)}.Data.Density;
        d2=Core{cn(2)}.Data.Density;
        if size(d1,2)>size(d1,1)
            d1 = d1';
        end
        if size(d2,2)>size(d2,1)
            d2 = d2';
        end
        indend=min(length(d1),length(d2));
        if indend == length(d1)
            depth = Core{cn(2)}.Data.Depth;
        else
            depth = Core{cn(1)}.Data.Depth;
        end
        if param.lag>=0
            densdiff = [zeros(param.lag, 1); d2(1:(indend-param.lag))-d1(1+param.lag:indend)];
            coef=corrcoef(d1(100+param.lag:indend),d2(100:(indend-param.lag)),'rows','pairwise');
        else
            densdiff = [d2(1-param.lag:indend)-d1(1:(indend+param.lag))];
            coef=corrcoef(d1(100:(indend+param.lag)),d2(100-param.lag:indend),'rows','pairwise');
        end
        r=coef(1,2);

        if strcmp(param.PlotDifference,'yes')  
            text=sprintf('lag = %i cm \n r = %f', -param.lag, r);
            h = annotation('textbox',...
            [0.05 0.1 0.25 0.08],...
            'Color','k',...
            'String',text,...
            'HorizontalAlignment','center',...
            'FontWeight','bold',...
            'FitBoxToText','off',...
            'EdgeColor','none');
  
            axes(ha(2));
            hold on
            plot((MySmooth(densdiff,param.span)), depth(1:length(densdiff)),'k','LineWidth',1.5)
            plot(((MySmooth(densdiff,param.span)+50)), depth(1:length(densdiff)),'Color',[0.5 0.5 0.5],'LineWidth',0.5)
            plot(((MySmooth(densdiff,param.span)-50)), depth(1:length(densdiff)),'Color',[0.5 0.5 0.5],'LineWidth',0.5)
            get(gca,'YLim')
            plot(zeros(2,1),get(gca,'YLim'),'--k')
            set(gca,'XMinorTick','on','YMinorTick','on',...
                'ydir','rev','YAxisLocation','right')
            axis tight
            ylim(ylimit)        
            legend('Density difference','Uncertainty',...
                'Location','NorthOutside')
            legend('boxoff')
            xlabel('Difference (kg/m^3)')
            box on
        end
    end
end