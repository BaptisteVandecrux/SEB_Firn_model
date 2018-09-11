function [minr, depth, stratmatch, bestlag] = StratStudy(Core, cn)

% transforming stratigraphy into numeric array 
%(firn = 0, ice_100% = 1, ice_x% =x/100)
    disp('Using StratStudy:')

    endind=min( length(Core{cn(1)}.Data.Depth), length(Core{cn(2)}.Data.Depth));
    stratmatch=[];
    numstrat(:,1)=NaN(endind,1);
    numstrat(:,2)=NaN(endind,1);

    for i=1:endind
        if strcmp(Core{cn(1)}.Data.Type{i},'ice') 
            numstrat(i,1)=Core{cn(1)}.Data.Type_perc(i);
        else
            numstrat(i,1) = 0;
        end
        
        if strcmp(Core{cn(2)}.Data.Type{i},'ice')                              
            numstrat(i,2)=Core{cn(2)}.Data.Type_perc(i);
        else
            numstrat(i,2) = 0;
        end
    end
        
    trylag=-50:50;
    
    r=[];
    for lag=trylag
        if lag>=0
            stratmatch = [zeros(lag, 1); numstrat(1+lag:endind,1)-numstrat(1:(endind-lag),2)];
        else
            stratmatch = numstrat(1:(endind+lag),1)-numstrat(1-lag:endind,2);
        end
        r =  [r sumabs(stratmatch(100:end))];
        depth=Core{cn(1)}.Data.Depth(1:endind);
    end
    [minr, ind] =min(r);
    bestlag=trylag(ind);

    figure 
    plot(trylag,r)
    title(sprintf('lag plot based on stratigraphy\n for %s and %s',Core{cn(1)}.Info.Name,Core{cn(2)}.Info.Name));

    fprintf('Best lag: %i cm (r = %f)\n',bestlag,minr);

    [~, ~, ~] = OverlapStrat(Core,cn,bestlag);
%     print(sprintf('%s - %s.png',Core{cn(j,1)}.Info.Name,Core{cn(j,2)}.Info.Name),'-dpng')
end