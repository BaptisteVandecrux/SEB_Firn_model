function [data_site] = CompactionAnalysis(DataCompaction, MetadataCompaction, ID)
   if ischar(ID)
       site = ID;
        site_num = find(strcmp(MetadataCompaction.sitename,site));

        disp('Instruments ID at that site:')
        instrument_list = MetadataCompaction.FC_instrument_IDs(site_num,MetadataCompaction.FC_instrument_IDs(site_num,:)~=0);
        disp(instrument_list)
        ID = instrument_list;
        data_site = CompactionAnalysis(DataCompaction,MetadataCompaction,ID);
   else
       for i = 1:length(ID)
           data_site{i} = DataCompaction{ID(i)};
            i_site=[];
            for ii =1:length(MetadataCompaction.sitename)
                temp = strfind(MetadataCompaction.sitename(ii),unique(data_site{i}.sitename));
                if ~isempty(temp{1})
                    i_site = ii;
                    break
                end
            end
            time_obs = datetime(datestr(data_site{i}.time));
            save{i}=time_obs;
            [ii, jj] = find(MetadataCompaction.FC_instrument_IDs==ID(i));
            
            ind_nan = isnan(data_site{i}.Compaction_Distance_m);

            % The "Compaction_Distance_m" given in Mike's files is the
            % length of a 2 m cable linked to the rewinding system. This
            % cable i further extended to a secondary cable down to the
            % bottom of the borehole. Therefor the borehole length reads as:            
            
            % Borehole Distance = Wire Measurement + (*Original* Borehole
            % Distance [metadata] - Original Wire Measurement at Day-0)
        
            data_site{i}.Compaction_Distance_m(~ind_nan) = ...
                data_site{i}.Compaction_Distance_m(~ind_nan) ...
                + (abs(MetadataCompaction.FC_borehole_initial_length_m(ii,jj)) ...
                - data_site{i}.Compaction_Distance_m(find(~ind_nan,1,'first')));
            data_site{i}.Compaction_Distance_m_1 = smooth(data_site{i}.Compaction_Distance_m,2*7,'lowess');
            data_site{i}.Compaction_Distance_m_2 = smooth(data_site{i}.Compaction_Distance_m,4*7,'lowess');

            data_site{i}.Compaction_Distance_m_1(ind_nan) = NaN;
            data_site{i}.Compaction_Distance_m_2(ind_nan) = NaN;
            a = data_site{i}.time(1:end-1)-data_site{i}.time(2:end);
            if sum(a>=0)>0
                fprintf('Duplicate time detected for instrument %i',ID(i))
                ind_err = find(a==0);
                format longg
                tab_err = table(data_site{i}.daynumber_YYYYMMDD(ind_err), ...
                    data_site{i}.daynumber_YYYYMMDD(ind_err +1), ...
                    data_site{i}.daynumber_YYYYMMDD(ind_err -1), ...
                    data_site{i}.Compaction_Distance_m(ind_err), ...
                    data_site{i}.Compaction_Distance_m(ind_err+1),...
                    'VariableName',{'Time_1','Time_2','Time_prev','Val_1','Val_2'})
            end
            
            data_site{i}.Compaction_Rate_md = ...
                [data_site{i}.Compaction_Distance_m(1:end-1)-data_site{i}.Compaction_Distance_m(2:end); 0];
            data_site{i}.Compaction_Rate_md_1 = ...
                [data_site{i}.Compaction_Distance_m_1(1:end-1)-data_site{i}.Compaction_Distance_m_1(2:end); 0];
            data_site{i}.Compaction_Rate_md_2 = ...
                [data_site{i}.Compaction_Distance_m_2(1:end-1)-data_site{i}.Compaction_Distance_m_2(2:end); 0];
            ind_err = ([data_site{i}.time(1:end-1)-data_site{i}.time(2:end); 1] ~= -1);
            data_site{i}.Compaction_Rate_md_2(ind_err) = NaN;
            data_site{i}.Compaction_Rate_md_1(ind_err) = NaN;
            data_site{i}.Compaction_Rate_md_0(ind_err) = NaN;
            
       end
       
   i_remove = [];
    for i = 1:length(ID)
        [ii, jj] = find(MetadataCompaction.FC_instrument_IDs==ID(i));
        time_obs = data_site{i}.time;

        if sum(~isnan(data_site{i}.Compaction_Distance_m_2))<10
            fprintf('%s - Instr. #%i - No data available',...
                MetadataCompaction.sitename{ii} , ID(i));  
            i_remove = [i_remove i];
        end
    end
    ID(i_remove)= [];
        
    f= figure('OuterPosition',[0 0 20 20*(length(ID)+2)/3]);
    ha = tight_subplot(length(ID)+2,1,0.07, 0.07, [0.12 0.08]);
  
    set(f,'CurrentAxes',ha(1))
    plot(data_site{1}.time,data_site{1}.AirTemp_median_C,'k','LineWidth',2)
    axis tight
    ylabel('Air temp-\newlineerature (^oC)','Interpreter','tex')
    set(gca,'XTickLabel',[],'XLim',[min(data_site{1}.time),max(data_site{1}.time)])
    set_monthly_tick(time_obs,gca);
    h_tit = title(MetadataCompaction.sitename{ii});
    
    set(f,'CurrentAxes',ha(2))
    ind_nonan = ~isnan(data_site{1}.SonicRangeDist_Corrected_m);
    ind=find(ind_nonan,1,'first');
    plot(data_site{1}.time,...
        hampel(data_site{1}.SonicRangeDist_Corrected_m(ind)-data_site{1}.SonicRangeDist_Corrected_m,14),...
        'k','LineWidth',2)
    axis tight
    ylabel('Surface \newlineHeight (m)','Interpreter','tex')
    set(gca,'XLim',[min(data_site{1}.time),max(data_site{1}.time)])
    set_monthly_tick(time_obs,gca);


    for i = 1:length(ID)
        
        [ii, jj] = find(MetadataCompaction.FC_instrument_IDs==ID(i));
        time_obs = data_site{i}.time;

        text_title = sprintf('Instr. #%i - Length: %0.2f m - Depth of top: %0.2f m',...
             ID(i),...
            -MetadataCompaction.FC_borehole_initial_length_m(ii,jj),...
            -MetadataCompaction.FC_borehole_top_from_surface_m(ii,jj));
        
        
    set(f,'CurrentAxes',ha(2+i))
            hold on
            [ax,h1, h2] = plotyy(data_site{i}.time,data_site{i}.Compaction_Distance_m_2,...
                data_site{i}.time,data_site{i}.Compaction_Rate_md_2.*1000);
            ylim(ax(2),[-0.5, 2])
            h1.LineWidth = 2;
            h2.LineWidth = 2;

            ind_nan = isnan(data_site{i}.Compaction_Distance_m);
            ind=find(~ind_nan,1,'first');

            axes(ax(1))
            hold on
            plot(data_site{i}.time(ind:(ind+60)), ...
                data_site{i}.Compaction_Distance_m_2(ind:(ind+60)), ...
                'Color',RGB('light light blue'),'LineWidth',2);
            set(gca,'XTickLabel',[],...
                'YLim',[min(data_site{i}.Compaction_Distance_m_2), max(data_site{i}.Compaction_Distance_m_2)],...
                'XLim',[min(data_site{1}.time),max(data_site{1}.time)])
            set_monthly_tick(time_obs,gca);

            axes(ax(2))
            hold on
            plot(data_site{i}.time, data_site{i}.Compaction_Rate_md.*1000, ...
                'Color',[0.8 0.8 0.8]);
            plot(data_site{i}.time,  data_site{i}.Compaction_Rate_md_2.*1000, ...
                'r','LineWidth',2);
            plot(data_site{i}.time(ind:(ind+60)), ...
                data_site{i}.Compaction_Rate_md_2(ind:(ind+60)).*1000, ...
                'Color',RGB('light red'),'LineWidth',2);

            set(gca,'XTickLabel',[],...
                'YLim',[-0.5 2],...
                'YTick',-0.5:0.5:2,...
                'YTickLabel',-0.5:0.5:2,...
                'XLim',[min(data_site{1}.time),max(data_site{1}.time)])
            set_monthly_tick(time_obs,gca);

            title(text_title)
            ax(1).YTickMode = 'auto';
            for kk = 1:2
                if size(ax(kk).YTickLabel,1)>3
                for k = 1:2:size(ax(kk).YTickLabel,1)
                    try ax(kk).YTickLabel{k} = '';
                    catch me
                        ax(kk).YTickLabel(k,:)=repmat(' ',1,size(ax(kk).YTickLabel,2));
                    end
                    
                end
                ax(kk).Box = 'off';
            end
            ax(2).XAxisLocation='top';
            ax(2).XTickLabel = '';
            ax(2).XAxis.Visible = 'on';
            if i == floor(length(ID)/2)
                h_1 = ylabel(ax(1),'Borehole length (m)','Interpreter','tex');
                h_1.Units = 'normalized';
                h_1.Position(1) = -0.1;
                h_1.Position(2) = -length(ID)/2;

                h_2 = ylabel(ax(2),'Compaction rate (mm/day)','Interpreter','tex');
                h_2.Units = 'normalized';
                h_2.Position(1) = 1.04;
                h_2.Position(2) = -length(ID)/2;
                
            end
        end
    end

    NameFile = sprintf('./Output/%s',...
        MetadataCompaction.sitename{ii});
    print(f,NameFile,'-dpng'); 
end