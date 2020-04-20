function [var_new] = ConvertToGivenDepthScale(depth_old, var_old, depth_new,opt)
% Interpolates the depth profile of a given variable old_var 
% (temperature grain size....) according to a given scale new_depth in real m

    transpose_at_the_end = 0;

    if isrow(depth_old) ~= isrow(depth_new)
        error('Old and new scale should be both rows or both columns.')
    elseif ~isrow(depth_old)
        transpose_at_the_end = 1;
       depth_old=depth_old';
       var_old=var_old';
       depth_new = depth_new';
    end
    var_new = NaN(size(depth_new));

    switch opt
        case 'linear'
        % the variable is interpolated on each stamp of the new scale
        var_new = interp1(depth_old,var_old,depth_new,'linear','extrap');
        
        case 'intensive'
            % intensive variables do not depend on the size of the system
            % f.e. if you know the density of a core section and you cut it
            % in two, the two sub-sections can be assigned the same density
            % as the original section. However we need to take into account
            % into the re-sampling the case when a new section is composed
            % of two sections in the old scale. Then the density of that
            % new section is the thickness-weighted average of the two
            % original sections.
            % example: depth_old = 1:4; density_old = 100:100:400;
            %           depth_new = [0.1 0.2 0.6 1.2 3.5];  
            % density_new = [100.0000  100.0000  100.0000  133.3333
            % 286.9565];

            if depth_new(end)>depth_old(end)
                depth_old = [depth_old, depth_new(end)];
                var_old = [var_old, var_old(end)];
            end
            left_neighbour_in_new = depth_new;
            left_neighbour_in_new(1) = 0;
            left_neighbour_in_new(2:end) = depth_new(1:end-1);

            left_neighbour_in_old =  interp1([0 depth_old],[0 depth_old],depth_new,'previous');
            
            ind_type_1 = left_neighbour_in_new >= left_neighbour_in_old;
            ind_type_2 = find(left_neighbour_in_new < left_neighbour_in_old);
            
            var_new(ind_type_1) = interp1([0 depth_old],...
                [var_old(1) var_old],...
                depth_new(ind_type_1),'next');
            
            depth_merged = [depth_old depth_new];
            var_merged = [var_old var_new];
            
            [depth_merged, ind_sorted] = sort(depth_merged);
            var_merged = var_merged(ind_sorted);
            var_merged(isnan(var_merged)) = interp1([0  depth_merged(~isnan(var_merged))],...
                [var_merged(1) var_merged(~isnan(var_merged))],...
                depth_merged(isnan(var_merged)),'next');

            thick_merged = depth_merged;
            thick_merged(2:end) = depth_merged(2:end) - depth_merged(1:end-1);
            
            for i = ind_type_2
                i_in_merged = discretize(depth_merged, ...
                    [left_neighbour_in_new(i) depth_new(i)]);
                i_in_merged(isnan(i_in_merged))= 0;
                if i~=1
                    i_in_merged(find(i_in_merged,1,'first')) = 0; % transforming the first one into 0
                end
                var_new(i) = sum(var_merged(i_in_merged==1).*thick_merged(i_in_merged==1)) ...
                    ./sum(thick_merged(i_in_merged==1));
            end
            
        case 'extensive'
            % extensive values depend on the size of the system
            % for example when downscaling the liquid water content of one
            % cell into two equally sized cells then the lwc in the new
            % cells are half of the lwc of the original ones
            % example: depth_old = 1:4; lwc_old = [1 0 0 1];
            %           depth_new = [0.1 0.2 0.6 1.2 3.5];  
            
            if depth_new(end)>depth_old(end)
                thick_last_old = depth_old(end)-depth_old(end-1);
                depth_old = [depth_old, depth_new(end)];
                thick_last_old_new = depth_old(end)-depth_old(end-1);
                var_old = [var_old, var_old(end)/thick_last_old*thick_last_old_new];
            end
            
            depth_merged = sort([depth_old depth_new]);
                        
            thick_merged = depth_merged;
            thick_merged(2:end) = depth_merged(2:end) - depth_merged(1:end-1);
            
            thick_old = depth_old;
            thick_old(2:end) = depth_old(2:end) - depth_old(1:end-1);

            ind_bin = discretize(depth_merged,[0 depth_old],'IncludedEdge','right');

            if sum(isnan(ind_bin))>1
                error('Some depths asked in depth_new not covered by depth_old')
            end
            var_merged = var_old(ind_bin).* thick_merged ./ thick_old(ind_bin);
            
            ind_bin_new = discretize(depth_merged,[0 depth_new],'IncludedEdge','right');
            var_merged(isnan(ind_bin_new)) =  [];
            ind_bin_new(isnan(ind_bin_new)) =  [];
            
            var_new = accumarray(ind_bin_new', var_merged')';
    end
    
    if transpose_at_the_end==1
       var_new=var_new';
    end
end