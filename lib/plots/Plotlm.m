function [lm, ah] = Plotlm(x,y,varargin)
param.Annotation = 'on';
param.PlotCF = 'off';
param.Unit = '^oC/dec';
param.Coeff = 1;
param.Color = 'k';
param.LineStyle = '--';
param.LineWidth = 1.5;
for i =1:2:length(varargin)
    param.(varargin{i}) = varargin{i+1};
end


lm = fitlm(x,y);
plot(x, ...
    lm.Coefficients.Estimate(2)*x+lm.Coefficients.Estimate(1),...
    'Color',param.Color,'LineWidth',param.LineWidth,...
    'LineStyle',param.LineStyle)
if strcmp(param.PlotCF , 'on')
    [~,ci1] = predict(lm,unique(x)',...
        'Alpha',0.05,'Simultaneous',true);
    plot(unique(x),ci1,':k');
end

if strcmp(param.Annotation,'on')
    slope = lm.Coefficients.Estimate(2)*param.Coeff;

    ah = annotation('textbox',...
        'String',sprintf('slope: %0.2f %s \n p-value: %0.2f',...
        slope*10,...
        param.Unit, ...
        max(lm.Coefficients.pValue)),...
        'LineStyle','-',...
        'FontSize',13);
    set(ah,'Parent',gca)
    set(ah,'position',[x(1)-3 ...
        lm.Coefficients.Estimate(2)*x(1)+lm.Coefficients.Estimate(1)-1 ...
        0.3 0.3])
else
    ah =0;
end
end
