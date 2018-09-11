function setmonth(source,event)
    Tsurf_obs = source.UserData.Tsurf_obs;
    Tsurf = source.UserData.Tsurf;
    month_obs = source.UserData.month_obs;
    c = source.UserData.c;
    
    if source.Value == 1
        month = month_obs;
    else
        month = source.Value -1;
    end

    scatter(Tsurf_obs(month_obs==month)-c.T_0,...
        Tsurf(month_obs==month),[],13 - month_obs(month_obs==month),'.')
    hold on
    plot(get(gca,'ylim'), get(gca,'ylim'),'k')

end