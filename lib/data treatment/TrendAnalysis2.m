function [trend, p_value, f] = TrendAnalysis2(periods, SEB_hour,varname, ...
    station, station2,opt1,opt2,vis,unit)
    trend = table;
    trend.no = (1:size(periods,1))';

    for i = 1:size(periods,1)
        trend.start_year(i) = periods(i,1);
        trend.end_year(i) = periods(i,2);
    end
    ymin=NaN; ymax = NaN;
    p_value=trend;

    f = figure('Visible',vis);
    ha = tight_subplot(3,3,[0.045 0.01],[0.2 0.08],0.07);
    for ii =1:length(station)   
        
        Time = SEB_hour{ii}.time;
        if strcmp(unit,'MJ')
            var = SEB_hour{ii}.(varname)/1000000;
        elseif strcmp(unit,'kJ')
            var = SEB_hour{ii}.(varname)/1000;
        else
            var = SEB_hour{ii}.(varname);
        end
        data = table(Time, var,'VariableNames',{'time','Var'});
        switch opt1
            case 'yearly'
                data_yearly = AvgTable(data,'yearly',opt2);
            case 'JJA'
                data_yearly = AvgTableJJA(data,opt2);
        end
        DV = datevec(data_yearly.time);
        years = DV(:,1);
        var_year = data_yearly.Var;

        set(f,'CurrentAxes',ha(ii))
        hold on
        scatter(years,var_year,'fill')
        [lm, ah] = Plotlm(years,var_year,'Annotation','off');
        h=plot(NaN,NaN,'w');
        lgd =legend(h,sprintf('slope: %0.2f %s \n p-value: %0.2f',...
            lm.Coefficients.Estimate(2)*10,...
            [unit '/dec'], ...
            coefTest(lm)),'Location','Best');
        legend boxoff
        lgd.FontSize = 12;

        title(station{ii})
        box on
        if ismember(ii,1:6)
           set(gca,'XTickLabel','')
       elseif ii==8
           xlabel('Year')
       end
       if ismember(ii,[2 3 5 6 8 9])
           set(gca,'YTickLabel','')
       elseif ii==4
            h_label = ylabel([varname,' ', opt1, ' ',opt2, ' (',unit,')'],'Interpreter','none');
       end
       xlim([1995 2017])
       set (gca,'XTick',1995:2014,'XMinorTick','off','YMinorTick','on','XTickLabelRotation',45)
        labels=get(gca,'XTickLabel');
        if ~isempty(labels)
        labels(1:2:end)={' '};
        set(gca,'XTickLabel',labels);
        end
        
        y_limit = get(gca,'YLim');
        ymax = max(ymax,y_limit(2));
        ymin = min(ymin,y_limit(1));
        
        name_slope = sprintf('slope_%s',station2{ii});
        name_pvalue = sprintf('pvalue_%s',station2{ii});
        trend.(name_slope) = (1:length(trend.no))';
        p_value.(name_pvalue) = (1:length(trend.no))';

        years_nonan = years(~isnan(var_year)); %years for which annual average is available

        for j = 1:length(trend.no)
            % we calculate things only if start and end years of the period are
            % available
            if ~isnan(trend.start_year(j))
                if ismember(trend.start_year(j),years_nonan) && ...
                        ismember(trend.end_year(j),years_nonan)
                    years_period = (trend.start_year(j):trend.end_year(j))';
                    var_period = var_year(ismember(years,years_period));

                    lm = fitlm(years_period, var_period);
                    trend.(name_slope)(j) = lm.Coefficients.Estimate(2)*10;
                    p_value.(name_pvalue)(j) = max(lm.Coefficients.pValue);
                else
                    trend.(name_slope)(j) = NaN;
                    p_value.(name_pvalue)(j) = NaN;
                end
            else
                years_period = years_nonan';
                var_period = var_year(ismember(years,years_period));

                lm = fitlm(years_period, var_period);
                trend.(name_slope)(j) = lm.Coefficients.Estimate(2)*10;
                p_value.(name_pvalue)(j) = coefTest(lm);
            end
        end
    end
    
    for ii =1:length(station)   
        set(f,'CurrentAxes',ha(ii))
        ylim([ymin ymax])
    end

    trend = table2array(trend);
    p_value = table2array(p_value);
end