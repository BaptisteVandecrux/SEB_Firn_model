function [r, depth, stratmatch] = OverlapStrat(Core, cn, lag) 
     set(0,'DefaulttextInterpreter','none')

    if length(cn)~=2
        disp('Enter only 2 indexes when using ''overlap'' option.')
        return
    end

    figure('units','normalized','outerposition',[0 0 0.4 1])
    for i=1:2
        if cn(i)<length(Core)
            disp(sprintf('Plotting core %i',cn(i)))

            if i == 1
                type=Core{cn(i)}.Data.Type;
                typeperc=Core{cn(i)}.Data.Type_perc;
                depth=Core{cn(i)}.Data.Depth;
            else
                type=Core{cn(i)}.Data.Type(1:(end));
                typeperc=Core{cn(i)}.Data.Type_perc(1:(end));
                if lag>=0
                    depth=Core{cn(i)}.Data.Depth((1+lag):end);
                else
                    depth=[(lag:0)'; Core{cn(i)}.Data.Depth];
                end
            end

            if i==1
                color = [0.7 0.8 1];
                edge1=1;
                edge2=475;
                linecolor = 'g';
            else
                color = [0.85 0.8 0.95];
                edge1=475;
                edge2=950;
                linecolor='b';
            end
            
            hold on

            for j=1:min(length(depth),length(type))-1
                if (strcmp(type{j},'ice')||strcmp(type{j},'ice layer')) & ~isnan(depth(j)) & ~isnan(depth(j+1))                              
                    if ~isnan(typeperc(j)) & isnumeric(typeperc(j))
                        rectangle('Position',[edge1 depth(j) edge2 (depth(j+1)-depth(j))],...
                            'FaceColor',color,'EdgeColor','none')
                    else
                        rectangle('Position',[edge1 depth(j) edge2 (depth(j+1)-depth(j))],...
                            'FaceColor',color,'EdgeColor','none')
                    end
                end
            end
        else
            disp('Index exceeding dataset size')
        end 
    
    if i==1
        color = [0.7 0.85 1];
        edge1=0.01;
        edge2=475;
        linecolor = 'b';
        posbox=[0.2 0.87 0.3 0.08];
    else
        color = [0.85 0.8 0.95];
        edge1=475;
        edge2=950;
        linecolor=[0.6 0 0.6];
        posbox=[0.5 0.87 0.3 0.08];
    end
            
        aux2=depth(~isnan(depth));
        rectangle('Position',[edge1 aux2(end) edge2 max((2000-aux2(end)),0)],...
            'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')

        set(gca,'ydir','rev')
        set(gca,'Position',[0.2 0.1 .6 .75])
        plottitle= sprintf('Name: %s \n Location: %s \n Elevation: %i m',...
            Core{cn(i)}.Info.Name,Core{cn(i)}.Info.NearestCodeLocation,Core{cn(i)}.Info.Coord1(3));
        % Create textbox
        h = annotation('textbox',...
            posbox,...
            'Color',linecolor,...
            'String',plottitle,...
            'HorizontalAlignment','center',...
            'FontWeight','bold',...
            'FitBoxToText','off',...
            'EdgeColor','none');
        set(h,'interpreter','none') % The underscore will no longer work as subscript
        xlabel('Density (kg/m^3)')
        ylabel('Depth (cm)')
        ylim([min(lag,0) 800])        
        xlim([0 950])
    end

    endind=min(length(Core{cn(1)}.Data.Depth), length(Core{cn(2)}.Data.Depth));
    stratmatch=[];
    temp1=[];
    temp2=[];
    
    for i=1:endind
        if (strcmp(Core{cn(1)}.Data.Type{i},'ice')||strcmp(Core{cn(1)}.Data.Type{i},'ice layer'))
            temp1=[temp1 1];
        else
            temp1=[temp1 0];
        end
        if(strcmp(Core{cn(2)}.Data.Type{i},'ice')||strcmp(Core{cn(2)}.Data.Type{i},'ice layer'))
            temp2=[temp2 1];
        else
            temp2=[temp2 0];
        end
    end
    if lag>=0
        stratmatch = [zeros(lag, 1)' temp1(1+lag:endind)-temp2(1:(endind-lag))];
        r =  sumabs(stratmatch);
    else
        stratmatch = [temp1(1:(endind+lag))-temp2(1-lag:endind)];
        r =  sumabs(stratmatch);
    end
    
    stairs(stratmatch.*250+250,Core{cn(1)}.Data.Depth(1:length(stratmatch))+1, 'LineWidth',1.5)
    depth=Core{cn(1)}.Data.Depth(1:endind);
end