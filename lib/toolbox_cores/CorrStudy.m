function [maxr, depth, densdiff, bestlag] = CorrStudy(Core, cn) 
    % Calculates the coefficient of correlation (r) and 
    % difference in density profiles (densdiff) between 
    % two cores (of indexes cn) with the second one shifted
    % by lag (in cm). It also gives back the depth vector to 
    % be used with densdiff.
disp('Using CorrStudy:')

    trylag = -150:150;
    r=[];
    
    for lag = trylag
        if length(cn)~=2
            disp('Enter only 2 indexes')
            return
        end   

        d1=Core{cn(1)}.Data.Density;
        d2=Core{cn(2)}.Data.Density;
        if size(d1,2)>size(d1,1)
            d1 = d1';
        end
        if size(d2,2)>size(d2,1)
            d2 = d2';
        end
        depth=Core{cn(1)}.Data.Depth;
        indend=min(length(d1),length(d2));
        if lag>=0
            densdiff = [zeros(lag, 1); d2(1:(indend-lag))-d1(1+lag:indend)];
            coef=corrcoef(d1(100+lag:indend),d2(100:(indend-lag)),'rows','pairwise');
        else
            densdiff = d2(1-lag:indend) - d1(1:(indend+lag));
            coef=corrcoef(d1(100:(indend+lag)),d2(100-lag:indend),'rows','pairwise');
        end

        if sum(sum(isnan(coef)))>0 || sum(size(coef))==2
            coef = [0, 0];
        end
            r = [r, coef(1,2)];
    end

    [maxr, ind] =max(r);
    bestlag=trylag(ind);
%     figure 
%     plot(trylag,r)
%     title(sprintf('lag plot for %s and %s',Core{cn(1)}.Info.Name,Core{cn(2)}.Info.Name));
%     
    fprintf('Best lag: %i cm (r = %f)\n',bestlag,maxr);
    [~, ~, ~] = OverlapPlot(Core, cn, bestlag);
%     print(sprintf('%s - %s.png',Core{cn(1)}.Info.Name,Core{cn(2)}.Info.Name),'-dpng')
end