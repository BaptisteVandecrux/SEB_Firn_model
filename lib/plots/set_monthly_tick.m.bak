function [] = set_monthly_tick(varargin)
    time = varargin{1};
    if length(varargin)==2
        ax = varargin{2};
    else
        ax = get(gca); 
    end
    if isnumeric(time)
        aux = datevec(time);
        temp = struct();
        temp.Year = aux(:,1);
        temp.Month = aux(:,2);
        temp.Day = aux(:,3);
    end
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    if (time(end)-time(1))>365*30
        ax.XTick = [datenum(temp.Year(1):10:temp.Year(end),1,1)];
        ax.XTickLabel = temp.Year(1):10:temp.Year(end);
        box on
        ax.XAxis.MinorTickValues = ...
            datenum(temp.Year(1):temp.Year(end));
        set(gca,'XMinorTick', 'on',...
            'YMinorTick','on',...
            'XTick', [datenum(temp.Year(1):10:temp.Year(end),1,1)],...
            'XTickLabel', temp.Year(1):10:temp.Year(end)); %,...
    %         'MinorTickValues', [datenum(temp.Year(1),temp.Month(1) +1:temp.Month(1)+12*(temp.Year(end)-temp.Year(1)+1),1)]);
        datetick('x','yyyy','keeplimits','keepticks')
%         if (time(end)-time(1))>365*8
%             xticklabels = get(gca,'XTickLabel');
%             for i =1:length(xticklabels)
%                 if floor(i/2)==i/2
%                     xticklabels(i,:) = '    ';
%                 end
%             end
%             set(gca,'XTickLabel',xticklabels); 
%             set(gca,'XTickLabelRotation',45);
%         end        
    elseif (time(end)-time(1))>365
        ax.XTick = [datenum(temp.Year(1):temp.Year(end),1,1)];
        ax.XTickLabel = temp.Year(1):temp.Year(end);
        box on
        ax.XAxis.MinorTickValues = ...
            [datenum(temp.Year(1),temp.Month(1) +1:temp.Month(1)+12*(temp.Year(end)-temp.Year(1)+1),1)];
        set(gca,'XMinorTick', 'on',...
            'YMinorTick','on',...
            'XTick', [datenum(temp.Year(1):temp.Year(end),1,1)],...
            'XTickLabel', temp.Year(1):temp.Year(end)); %,...
    %         'MinorTickValues', [datenum(temp.Year(1),temp.Month(1) +1:temp.Month(1)+12*(temp.Year(end)-temp.Year(1)+1),1)]);
        datetick('x','yyyy','keeplimits','keepticks')
        if (time(end)-time(1))>365*8
            xticklabels = get(gca,'XTickLabel');
            for i =1:length(xticklabels)
                if floor(i/2)==i/2
                    xticklabels(i,:) = '    ';
                end
            end
            set(gca,'XTickLabel',xticklabels); 
            set(gca,'XTickLabelRotation',45);
        end

    else
%         if temp.Month(end)>=temp.Month(1)
%             XTicks = [datenum(temp.Year(1),temp.Month(1):temp.Month(end),1)];
%         else
            XTicks = [datenum(temp.Year(1),temp.Month(1):12,1) datenum(temp.Year(1)+1,1:temp.Month(end),1)];
%         end

           set(gca,'XMinorTick', 'on',...
            'YMinorTick','on',...
            'XTick', XTicks); 
            datetick('x','dd-mmm-yyyy','keeplimits','keepticks')

    end
end