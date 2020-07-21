function [] = Plot_seasonal_temp(data,seas,col,FS)

    fc_name = sprintf('AvgTable%s',seas);
    data_seas= feval(fc_name, data,'nanmean');
    
    T_seas = data_seas.Var;
    stairs([data_seas.time; data_seas.time(end)+365],...
    [T_seas; T_seas(end)],'Color',col,'LineWidth',2)
    DV = datevec(data_seas.time);
    years = DV(:,1);
    lm = fitlm(data_seas.time,T_seas);

    switch seas
        case 'JJA'
            y_txt = 2003.8;
        case 'SON'
            y_txt = 2005.8;
        case 'DJF'
            y_txt = 2008;
        case 'MAM'
            y_txt = 2010.5;
    end
    if lm.Coefficients.Estimate(2)*365*10 > 0.1
        text_lm = sprintf('%+0.2f',...
            lm.Coefficients.Estimate(2)*365*10);
    else
        text_lm = sprintf('< 0.1',...
            lm.Coefficients.Estimate(2)*365*10);
    end
    text(datenum(y_txt,1,1),0,...
        text_lm,...
        'Color',col,'FontWeight','bold',...
        'Interpreter','tex','FontSize',FS)
    
end