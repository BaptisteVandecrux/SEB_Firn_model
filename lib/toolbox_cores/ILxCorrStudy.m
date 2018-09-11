function [r, depth, densdiff, bestlag] = ILxCorrStudy(Core, cn)
% deos the same as CorrStudy but using the xcorr function on scaled density
% profiles
disp('Using xCorrStudy:')
%scaling density profiles
    T1 = Core{cn(1)}.Data.Type_perc;
    scaled1 = (T1 - nanmean(T1)) / nanstd(T1);
    scaled1(isnan(T1)) = 0;

    T2 = Core{cn(2)}.Data.Type_perc;
    scaled2 = (T2 - nanmean(T2)) / nanstd(T2);
    scaled2(isnan(T2)) = 0;

    [C, lag] = xcorr(scaled1(100:end),scaled2(100:end));

%     figure
%     plot(lag,C,'k');
%     ylabel('Amplitude');
%     grid on
%     title(sprintf('Cross-correlation between %s and %s',...
%         Core{cn(1)}.Info.Name, Core{cn(2)}.Info.Name))
%     xlabel('lag(cm)');

    [maxr,I] = max(abs(C));
    bestlag = lag(I);

    fprintf('Best lag: %i cm (r = %f)\n',bestlag,maxr);
    [r, depth, densdiff] = OverlapPlot(Core, cn, bestlag);
end