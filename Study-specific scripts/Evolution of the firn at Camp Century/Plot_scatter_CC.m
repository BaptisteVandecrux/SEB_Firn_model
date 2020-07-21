function [f] = Plot_scatter_CC(data_surf,VarList, data_AWS_GITS,  ...
    data_AWS_CEN,VarList2, vis,LabelList, NameFile,model_list, word)

ABC = char(65:90) ;
%%
f = figure('Visible',vis,'OuterPosition',[0 0 30 40],'DefaultAxesFontSize',10);
ha = tight_subplot(4,5, [0.1 0.04], [0.07 0.13], 0.07);
count = 0;
for ii = 1:4
    
%     data_surf{ii}.time  = data_surf{ii}.time(1) + ...
%         ((1:length(data_surf{ii}.time))' - 1)/24;
    ind_GITS = find(data_AWS_GITS.time<data_surf{ii}.time(end)+0.001);
    ind_CEN = find(data_AWS_CEN.time<data_surf{ii}.time(end)+0.001);

    ind_mod_GITS = find(and(data_surf{ii}.time>data_AWS_GITS.time(1)-0.001,...
        data_surf{ii}.time<data_AWS_GITS.time(end)+0.001));
    ind_mod_CEN = find(and(data_surf{ii}.time>data_AWS_CEN.time(1)-0.001,...
        data_surf{ii}.time<data_AWS_CEN.time(end)+0.001));
    clearvars h

    col_aws = lines(2);

    fac =1;
    if exist('opt','var')==1
        if isnumeric(opt)
            fac=opt;
            clearvars opt
        end
    end
    for k = 1:length(VarList)
        count = count+1;
        set(f,'CurrentAxes',ha(count))
        hold on

        if ismember(VarList2{k},data_AWS_GITS.Properties.VariableNames)
% 
%             [data_AWS_GITS.(VarList2{k})(ind_GITS), ~] = ...
%                 temporal_shift(data_AWS_GITS.(VarList2{k})(ind_GITS),...
%                 data_surf{ii}.(VarList{k})(ind_mod_GITS)*fac,...
%                 VarList2{k},2);
            
            plot(data_AWS_GITS.(VarList2{k})(ind_GITS),...
                data_surf{ii}.(VarList{k})(ind_mod_GITS)*fac,...
                '.','MarkerEdgeColor',col_aws(1,:),...
                'MarkerFaceColor',col_aws(1,:));
        end
%         [data_AWS_CEN.(VarList2{k})(ind_CEN), ~] = ...
%             temporal_shift(data_AWS_CEN.(VarList2{k})(ind_CEN),...
%             data_surf{ii}.(VarList{k})(ind_mod_CEN)*fac,...
%             VarList2{k},2);
        plot(data_AWS_CEN.(VarList2{k})(ind_CEN),...
            data_surf{ii}.(VarList{k})(ind_mod_CEN)*fac,'.',...
            'MarkerEdgeColor',col_aws(2,:),...
            'MarkerFaceColor',col_aws(2,:));

        axis tight ; box on
        plot([min(ha(k).YLim(1),ha(k).XLim(1)), max(ha(k).YLim(2),ha(k).XLim(2))],...
            [min(ha(k).YLim(1),ha(k).XLim(1)), max(ha(k).YLim(2),ha(k).XLim(2))],...    
            'k')
        
        if ismember(VarList2{k},data_AWS_GITS.Properties.VariableNames)
            lm1 = fitlm(data_AWS_GITS.(VarList2{k})(ind_GITS),...
                data_surf{ii}.(VarList{k})(ind_mod_GITS)*fac);
            text1 = sprintf(['\\color[rgb]{0, 0.447, 0.741}'...
                'R^2 = %0.2f RMSE = %0.1f ME = %0.1f'],...
                lm1.Rsquared.Ordinary, lm1.RMSE,...
                nanmean( data_surf{ii}.(VarList{k})(ind_mod_GITS) - ...
                data_AWS_GITS.(VarList2{k})(ind_GITS)));
        else
            text1 = '\color[rgb]{0, 0.447, 0.741}Not available at GITS';
        end
      
        lm2 = fitlm(data_AWS_CEN.(VarList2{k})(ind_CEN),...
            data_surf{ii}.(VarList{k})(ind_mod_CEN)*fac);
        
        text2 = sprintf(['\\color[rgb]{0.8500, 0.3250, 0.0980}'...
            'R^2 = %0.2f RMSE = %0.2f ME = %0.2f'],...
            lm2.Rsquared.Ordinary, lm2.RMSE,...
            nanmean(data_surf{ii}.(VarList{k})(ind_mod_CEN) - ...
            data_AWS_CEN.(VarList2{k})(ind_CEN)));
        
        title_text = {LabelList{k};  text1; text2};
        h_tit = title(title_text,'Interpreter','tex',...
            'FontSize',9,...
            'HorizontalAlignment','center');
        h_tit.Units = 'normalized';
        h_tit.Position(2) = h_tit.Position(2)*1.05;
        h_now = gca;
        h_text = text(h_now.XLim(2)-(h_now.XLim(2)-h_now.XLim(1))*0.2,...
            h_now.YLim(1)+(h_now.YLim(2)-h_now.YLim(1))*0.15, ABC(count));
        h_text.FontWeight='bold';
        h_text.FontSize=13;
        
        if k == 1
            ylabel({[word,' values from '], model_list{ii}},'FontSize',9);
        end
    end

end
h_x = xlabel(ha(18),'AWS observations');

h(1) = plot(NaN,NaN,'o',...
    'MarkerEdgeColor',col_aws(1,:),...
    'MarkerFaceColor',col_aws(1,:));
h(2) = plot(NaN,NaN,'o',...
    'MarkerEdgeColor',col_aws(2,:),...
    'MarkerFaceColor',col_aws(2,:));
legendflex(h,{'GITS AWS','CEN AWS'}, ...
    'ref', gcf, ...
       'anchor',  {'n' 'n'} , ...
       'buffer',[0 -5], ...
       'nrow',1, ...
       'fontsize',10,'Interpreter','none');

orient(f,'portrait')
print(f,NameFile,'-djpeg'); 
%% Second plot
% for ii = 1:2
%     data_surf{ii}.time  = data_surf{ii}.time(1) + ...
%         ((1:length(data_surf{ii}.time))' - 1)/24;
%     ind_GITS = find(ismember(data_AWS_GITS.time,data_surf{ii}.time));
%     ind_CEN = find(ismember(data_AWS_CEN.time,data_surf{ii}.time));
%     ind_mod_GITS = find(ismember(data_surf{ii}.time,data_AWS_GITS.time));
%     ind_mod_CEN = find(ismember(data_surf{ii}.time,data_AWS_CEN.time));
%     clearvars h
% 
%     col_aws = lines(2);
%     h(1) = plot(NaN,NaN,'o','MarkerEdgeColor',col_aws(1,:),...
%         'MarkerFaceColor',col_aws(1,:));
%     h(2) = plot(NaN,NaN,'o','MarkerEdgeColor',col_aws(2,:),...
%         'MarkerFaceColor',col_aws(2,:));
% 
%     
%     f = figure('Visible',vis,'OuterPosition',[0 0 25 16]);
%     ha = tight_subplot(4,2,...
%         0.05, 0.1, 0.12);
%     fac =1;
%     if exist('opt','var')==1
%         if isnumeric(opt)
%             fac=opt;
%             clearvars opt
%         end
%     end
%     for k = 1:length(VarList)
%         set(f,'CurrentAxes',ha(k))
%        
%         ylabel(LabelList{k});
% 
%         h_now = gca;
%         h_text = text(h_now.XLim(2)-(h_now.XLim(2)-h_now.XLim(1))*0.2,...
%             h_now.YLim(1)+(h_now.YLim(2)-h_now.YLim(1))*0.15, ABC(k));
%         h_text.FontWeight='bold';
%         h_text.FontSize=13;
% 
%         hold on
%         if ismember(VarList2{k},data_AWS_GITS.Properties.VariableNames)
%             plot(data_AWS_GITS.time(ind_GITS),...
%                 data_AWS_GITS.(VarList2{k})(ind_GITS));
%         end
%     
%         plot(data_AWS_CEN.time(ind_CEN),data_AWS_CEN.(VarList2{k})(ind_CEN));
%         plot(data_surf{ii}.time((ind_mod_GITS(1):ind_mod_CEN(end))),...
%             data_surf{ii}.(VarList{k})(ind_mod_GITS(1):ind_mod_CEN(end)));
%         set_monthly_tick(data_surf{ii}.time(ind_mod_GITS(1):ind_mod_CEN(end)));
% 
%     end
% 
%     h_x = xlabel('AWS observation');
%     h_x.Units = 'normalized';
%     h_x.Position(1) = h_x.Position(1)-2;
%     h_y = ylabel(['Adjusted values from ' model_list{ii}]);
%     h_y.Units = 'normalized';
%     h_y.Position(1) = h_y.Position(1)-4.1;
%     h_y.Position(2) = h_y.Position(2)+1;
%     % legendflex(h,{'GITS AWS','CEN_AWS'}, ...
%     %     'ref', gcf, ...
%     %        'anchor',  [2 6] , ...
%     %        'buffer',[0 -30], ...
%     %        'nrow',1,'box','off', ...
%     %        'fontsize',13,'Interpreter','none');
% 
%     orient(f,'portrait')
%     print(f,[NameFile(1:end-7) '_overlap_' model_list{ii}],'-dpng'); 
% end

end