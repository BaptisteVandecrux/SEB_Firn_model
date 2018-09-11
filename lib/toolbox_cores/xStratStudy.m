function [minr, depth, stratmatch, bestlag] = xStratStudy(Core, cn)

% transforming stratigraphy into numeric array 
%(firn = 0, ice_100% = 1, ice_x% =x/100)

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
    
    disp('Using xStratStudy:')
    %scaling density profiles
    T1 = numstrat(:,1);
    scaled1 = (T1 - nanmean(T1)) / nanstd(T1);
    scaled1(isnan(T1)) = 0;

    T2 = numstrat(:,2);
    scaled2 = (T2 - nanmean(T2)) / nanstd(T2);
    scaled2(isnan(T2)) = 0;

    [C, lag] = xcorr(scaled1(100:end),scaled2(100:end));   
    
    figure
    plot(lag,C,'k');
    ylabel('Amplitude');
    grid on
    title(sprintf('Cross-correlation between %s and %s (strat)',...
        Core{cn(1)}.Info.Name, Core{cn(2)}.Info.Name))
    xlabel('lag(cm)');

    [~,I] = max(abs(C));
    bestlag = lag(I);
    [minr, depth, stratmatch] = OverlapStrat(Core, cn, bestlag);
    fprintf('Best lag (strat): %i cm (r = %f)\n',bestlag,minr);
%     print(sprintf('%s - %s.png',Core{cn(j,1)}.Info.Name,Core{cn(j,2)}.Info.Name),'-dpng')
end