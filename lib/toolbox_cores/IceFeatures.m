function [ice_feat, depth_ice] = IceFeatures(Core, cn)
    depth = Core{cn}.Data.Depth;
    type = Core{cn}.Data.Type;
    type_perc = Core{cn}.Data.Type_perc;
    if ~isempty(Core{cn}.Data.Comment)
        comment = Core{cn}.Data.Comment;
    else
        comment=num2cell(NaN(size(type)));
    end
    
    depth_ice=[];
    if isempty(type)
        fprintf('\n No stratigraphy available for core %s.\n',...
            Core{cn}.Info.Name);
        ice_feat=[];
    else
        count = 1; % count number for ice features
        count_ice = 0; % count ice feature thickness
        
        %initializing
        ice_struct{count}.start_depth=NaN;
        ice_struct{count}.end_depth=NaN;
        ice_struct{count}.thickness=NaN;
        
        for i=2:length(type_perc)
            if type_perc(i)>0
                depth_ice = [depth_ice, depth(i)];
                if (~isempty(strfind(comment{i}, 'lens')) ...
                        && isempty(strfind(comment{i-1}, 'lens')) ...
                        && isempty(strfind(comment{i+1}, 'lens')))
                % this is for isolate ice lenses
                % ie "comment" contains "lens" but cell abovee and below 
                % dont contain "lens" or if it is an isolated ice layer

                    %by default the type_perc is used to define the ice thickness
                    % (50% ice = 5mm ice lense)
                    thick = Core{cn}.Data.Type_perc(i)/100;            

                    if ~isempty(strfind(comment{i}, 'mm'))
                    %if the comment contains the sub-cm thickness of the layer then
                    %it replaces what set above
                  
                        ind = strfind(comment{i}, 'mm');
                        if str2double(comment{i}(ind-1))~=0
                            thick = str2double(comment{i}(ind-1))/10;
                        else
                            thick = 1;
                        end
                    end

                    ice_struct{count}.start_depth = depth(i);
                    ice_struct{count}.thickness = thick;
                    ice_struct{count}.end_depth = depth(i)+thick;

                    count = count +1; %next ice layer will be counted as different
                    count_ice = 0; %reset ice thickness
                else
                    if  type_perc(i) >= 50
                        %to discard pipes and only work with ice layers extending 
                        %horizontally, what is composed of less than 50% ice is not
                        %counted as ice layer.

                        if count_ice == 0
                            %if no ice counted in count_ice
                            %then it is the beginning of new ice structure
                            ice_struct{count}.start_depth = depth(i);
                        end

                        count_ice = count_ice + 1;
                        if i<length(type)-1
                            if type_perc(i+1) < 50
                                %if no ice after or ice lower than 50 %
                                % then end of ice structure
                                ice_struct{count}.end_depth = depth(i);
                                ice_struct{count}.thickness = count_ice;
                                count_ice = 0; %reset ice thickness
                                count = count + 1; %next ice feature is distinct
                            end
                        else
                            if count_ice == 0
                                %if no ice counted in count_ice
                                %then it is the beginning of new ice structure
                                ice_struct{count}.start_depth = depth(i)
                            end

                            count_ice = count_ice + 1;

                            % end of ice structure
                            ice_struct{count}.end_depth = depth(i)
                            ice_struct{count}.thickness = count_ice;
                        end
                    else
                        disp('An ice feature with <50% ice was ignored.')
                    end
                end
            end
        end

    % Reorganizing the struct into a table
        ice_feat = zeros(length(ice_struct), 3);
        for i = 1:length(ice_struct)
            ice_feat (i,1) = ice_struct{i}.start_depth;
            ice_feat (i,2) = ice_struct{i}.end_depth;
            ice_feat (i,3) = ice_struct{i}.thickness;
        end

        ice_feat = array2table(ice_feat, 'VariableNames', ...
            {'start_depth', 'end_depth', 'thickness'});
        fprintf('\n Ice features of core: %s\n', Core{cn}.Info.Name);
        disp(ice_feat)

        
        figure
         set(0,'DefaulttextInterpreter','none')
        subplot(1,2,1)
        h = histogram(ice_feat.thickness,0:5:300);
        h.FaceColor = [0.8 0.5 0.5];
        h.EdgeColor = 'k';
        xlabel('Thickness of ice feature (cm)')
        ylabel('Frequency')
        ylim([0 15])
        
        subplot(1,2,2)
        h2 = histogram(depth_ice,0:200:2000);
        h2.FaceColor = [0.5 0.5 0.8];
        h2.EdgeColor = 'k';
        xlabel('Depth of ice feature (cm)')
        ylabel('Frequency')
        ylim([0 150])
        hold on 
        r=rectangle('Position', [length(type)+10 0.2...
            1980-length(type) 149]);
        r.EdgeColor='none';
        r.FaceColor = [0.8 0.8 0.8];
        
        h = suptitle(sprintf(['%s - %s - %0.1f %% of core is ice.\n',...
        'Nbr. ice features: %i. Mean thick. of feat.: %0.1f cm.'],...
            Core{cn}.Info.NearestCodeLocation, Core{cn}.Info.Name, ...
            length(depth_ice)/length(type)*100,...
            size(ice_feat,1), mean(ice_feat.thickness)));
        h.FontSize = 12;
    end
end