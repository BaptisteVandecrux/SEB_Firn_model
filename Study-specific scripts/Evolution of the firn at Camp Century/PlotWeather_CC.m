function [f,ha] = PlotWeather_CC(data,vis,VarList,LabelList,...
    NameFile,model_list,opt)
ABC = char(65:90) ;
ftsize = 12;
f = figure('Visible',vis,...
    'OuterPosition',[0 0 24 20*(length(VarList)+1)/8],...
    'DefaultAxesFontSize',ftsize);

if length(VarList)>4
    ha = tight_subplot(length(VarList)/2,2,...
        0.025, [0.1 0.05], 0.12);
else
    ha = tight_subplot(length(VarList)/2,2,...
        0.025, [0.15 0.1], 0.1);
end
fac =1;
        if exist('opt','var')==1
if isnumeric(opt)
    fac=opt;
    clearvars opt
end
        end
for k = 1:length(VarList)
            set(f,'CurrentAxes',ha(k))
            hold on
    for ii = 1:4    
        if exist('opt','var')==1
            data{ii}.(VarList{k})(isnan(data{ii}.(VarList{k}))) = 0;
            data{ii}.(VarList{k}) = cumsum(data{ii}.(VarList{k}));
        end 
        switch VarList{k}
            case 'Tsurf'
                data{ii}.(VarList{k})=data{ii}.(VarList{k})-273.15;
            case 'melt_mweq'
                data{ii}.(VarList{k})=data{ii}.(VarList{k})*1000;
        end
        stairs(data{ii}.time, data{ii}.(VarList{k})*fac,'LineWidth',1.5)
        DV = datevec(data{ii}.time);
        yr =DV(:,1);
        lm1 = fitlm(yr(yr<2020), data{ii}.(VarList{k})(yr<2020)*fac);
        if sum(yr>=2020)>1
            lm2 = fitlm(yr(yr>=2020), data{ii}.(VarList{k})(yr>=2020)*fac);
        
            fprintf('%s\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n', ...
                VarList{k},...
                nanmean(data{ii}.(VarList{k})(yr<2020)*fac),...
                nanstd(data{ii}.(VarList{k})(yr<2020)*fac),...
                lm1.Coefficients.Estimate(2)*10,...
                max(lm1.Coefficients.pValue),...
                nanmean(data{ii}.(VarList{k})(yr>=2020)*fac),...
                nanstd(data{ii}.(VarList{k})(yr>=2020)*fac),...
                lm2.Coefficients.Estimate(2)*10,...
                max(lm2.Coefficients.pValue));
        else
            fprintf('%s\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t-\t-\t-\t-\n', ...
                VarList{k},...
                nanmean(data{ii}.(VarList{k})(yr<2020)*fac),...
                nanstd(data{ii}.(VarList{k})(yr<2020)*fac),...
                lm1.Coefficients.Estimate(2)*10,...
                max(lm1.Coefficients.pValue));
        end
    end
    set_monthly_tick(datenum(1966:1/365:2101,1,1));
    axis tight
    if k/2 ==floor(k/2)
        set(gca,'YAxisLocation','right')
    end
    ylabel(LabelList{k},'Interpreter','tex')
    if ~ismember(k, [length(VarList), length(VarList)-1])
        set(gca,'XTickLabel','')
    else
            xticklabels = get(gca,'XTickLabel');
            for i =1:length(xticklabels)
                if floor(i/2)~=i/2
                    xticklabels(i,:) = '    ';
                end
            end
            set(gca,'XTickLabel',xticklabels); 
%             set(gca,'XTickLabelRotation',45);
    end

    ylimit = get(gca,'Ylim');
    tmp1 = nanmax(data{1}.(VarList{k})(1:10)*fac) + (ylimit(2)-ylimit(1))*0.5;
    tmp4 = nanmax(data{4}.(VarList{k})(1:10)*fac) + (ylimit(2)-ylimit(1))*0.5;
    set(gca,'Ylim',[ylimit(1) max([ylimit(2),tmp1, tmp4])]);
    if strcmp(VarList{k},'albedo')
        ylim([0.83, 0.89])
    end
    
    h_now = gca;
    h_text = text(h_now.XLim(1)+(h_now.XLim(2)-h_now.XLim(1))*0.02,...
        h_now.YLim(2)-(h_now.YLim(2)-h_now.YLim(1))*0.15, ABC(k));
    h_text.FontWeight='bold';
    h_text.FontSize=13;
    
    ha(k).YAxis.Exponent = 0;
    grid on;
end
ha(end).XLabel.String='Year';
ha(end-1).XLabel.String='Year';
legendflex(model_list, 'ref', gcf, ...
       'anchor',  [2 6] , ...
       'buffer',[0 -30], ...
       'nrow',1,'box','off', ...
       'fontsize',ftsize,'Interpreter','none');
   
orient(f,'portrait')
print(f,NameFile,'-djpeg'); 
end