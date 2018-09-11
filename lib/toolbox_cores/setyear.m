function setyear(source,event)
    T_avg = source.UserData.T_avg;
    PV_10m = source.UserData.PV_10m;
    ind_bin = source.UserData.ind_bin;
    my_fit = {source.UserData.my_fit};
    
    if source.Value == 1
        bin = ind_bin;
    else
        bin = source.Value -1;
    end
col = lines(max(ind_bin));
    scatter(T_avg(ind_bin==bin),...
        PV_10m(ind_bin==bin),100,...
        col(ind_bin(ind_bin==bin),:),'fill')
    hold on
        xlimit=get(gca,'XLim');
        x = xlimit(1):xlimit(end);
    for i = bin
%         y = my_fit{i}(1).*x.^3 + my_fit{i}(2).*x.^2 + my_fit{i}(3).*x +my_fit{i}(4);

        y = slmeval(x,my_fit{i});
        if i==1 || i ==5
            y(x>-14) = NaN;
        end
        plot(x,y,'Color',col(i,:),'Linewidth',1.5)
    end
        ylim([0 6])        
        xlim([-28 -10])
end