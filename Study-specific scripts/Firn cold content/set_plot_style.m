function [h_text] = set_plot_style(ii, station, time_mod, ylab)

switch ii
    case 2 
        station{ii} = "Dye-2";    
    case 6 
        station{ii} = "South Dome";
    case 9
        station{ii} = "Tunu-N";
end
if ~iscell(time_mod)
    tmp = cell(1,length(station));
    tmp{ii} =time_mod;
    time_mod = tmp;
end

    axis tight
    box on 
    time_ext = datenum(1998,1,1):1/24:datenum(2017,12,31);
%     set_monthly_tick(time_ext);
    set(gca,'XTick',datenum(2000:5:2020,1,1));
datetick('x','keeplimits','keepticks')
    h_text = text(time_ext(24*20), ...
       -20,...
        sprintf('%s %s',char(ii+96),station{ii}));
    h_text.FontSize = 15;
    h_text.FontWeight = 'bold';
    h_text.Units = 'Normalized';
    
    h_text.Position(1:2) = [0.03 0.85];
    xlim([datenum(1998,1,1) datenum(2017,12,31)])

    if ismember(ii,1:6)
       set(gca,'XTickLabel','')
    else
       xlabel('Year')
%        set_monthly_tick(time_mod{ii})
    end

   xlim([datenum(1998,1,1) datenum(2017,12,31)])
   if ismember(ii,[2 3 5 6 8 9])
%        set(gca,'YTickLabel','')
   end
   if ii==4
        ylabel(ylab,'Interpreter','tex')
   end
   grid on
end