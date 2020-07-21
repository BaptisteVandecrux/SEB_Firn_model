function [] = PlotRemoval(data1, data3, var_name_uni, station, vis,txt)
ind_remove = [];
for i = 1:length(var_name_uni)
    if ~ismember(var_name_uni{i},data1.Properties.VariableNames)
        ind_remove = [ind_remove i];
    end
end
var_name_uni(ind_remove) = [];

f = figure('Visible',vis);
ha = tight_subplot(length(var_name_uni),1,0.02, [0.05 0.02],[0.07 0.25]);
for i = 1:length(var_name_uni)
        set(f,'CurrentAxes',ha(i))
        hold on
        plot(data1.time,data1.(var_name_uni{i}),'Color',RGB('red'))
        plot(data3.time,data3.(var_name_uni{i}),'Color',RGB('dark green'))
    if i == 1
        h_tit = title(station);
        h_tit.Units = 'normalized';
        h_tit.Position(2) = h_tit.Position(2)-0.5;
    end
            xlabel('')
            ylabel('')

            if i~=length(var_name_uni)
                set(gca,'XTickLabel','');
            end
            h_leg = legend('Erroneous data',...
                sprintf('%s',var_name_uni{i}),...
                'Location','NorthWest');
            h_leg.Position(1) = h_leg.Position(1) +  0.69;

            set_monthly_tick(data1.time)
            if i <length(var_name_uni)
                set(gca,'XTickLabel','')
            end
                        axis tight
        end
        set(gca,'XTickLabelRotation',0)
        print(f,['./Output/Corrected/Plots/ErroneousData_' station txt],'-dpng')
