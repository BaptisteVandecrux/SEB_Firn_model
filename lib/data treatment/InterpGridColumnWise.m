function [depth_alt, M_alt, avg_shallow] = ...
    InterpGridColumnWise(depth_grid, M, d, opt, opt2, rho_alt)
    depth_alt = [depth_grid; d*ones(1,size(depth_grid,2))];

    % sorting so that NaN are located at depth 20 m
    [depth_alt, ind_sorted] = sort(depth_alt,1);
    
    M_alt = [M; NaN(1,size(M,2))];
    if size(M,2)==1
        M_alt = M_alt(ind_sorted);
    else
        for i = 1:size(M,2)
            M_alt(:,i) = M_alt(ind_sorted(:,i),i);
        end
    end
    ind_fixed_depth = isnan(M_alt);
    ind_all_nan =sum(ind_fixed_depth,1)>1;
    ind_fixed_depth(:,ind_all_nan) = 0;
    ind_fixed_depth(end,ind_all_nan) = 1;
    % we then interpolate at each time step
    ind_next =  boolean([zeros(1,size(M_alt,2)); ind_fixed_depth(1:end-1,:)]);
    ind_prev =  boolean([ind_fixed_depth(2:end,:); zeros(1,size(M_alt,2))]);
    
    % index of time steps where deepest thermistor is closest to fixed
    % depth
    ind_extrap = ind_fixed_depth(end,:) == 1;
    % in these case we use the last two elements to extrapolate to demanded
    % depth
    ind_next(:,ind_extrap)=0;
    ind_next(end-1,ind_extrap)=1;
    ind_prev(:,ind_extrap)=0;
    ind_prev(end-2,ind_extrap)=1;
    %same thing when we need to extrapolate towards the surface
    ind_extrap = ind_fixed_depth(1,:) == 1;
    % in these case we use the first two elements to extrapolate to demanded
    % depth
    ind_next(:,ind_extrap)=0;
    ind_next(1,ind_extrap)=1;
    ind_prev(:,ind_extrap)=0;
    ind_prev(2,ind_extrap)=1;
    
    switch opt
        case 'next'
            M_alt(ind_fixed_depth) = M_alt(ind_next);
        case 'lin'
            x1 = depth_alt(ind_prev);
            x2 = depth_alt(ind_next);
            y1 = M_alt(ind_prev);
            y2 = M_alt(ind_next);
            M_alt(ind_fixed_depth) = (y2-y1)./(x2-x1) .* (d -x1) + y1;
    end
    
    thickness_alt  = depth_alt;
    thickness_alt(2:end,:)  = depth_alt(2:end,:)-depth_alt(1:end-1,:);

    ind_below_thresh = depth_alt>d;
    M_shallow = M_alt;
    M_shallow(ind_below_thresh) = NaN;
    thickness_shallow = thickness_alt;
    thickness_shallow(ind_below_thresh) = NaN;
    
    switch opt2
        case 'MassWeighted'
            mass_shallow = rho_alt.*thickness_shallow;
            avg_shallow = nansum(M_shallow.*mass_shallow,1)...
                ./nansum(mass_shallow,1);
        case 'VolumeWeighted'
            avg_shallow = nansum(M_shallow.*thickness_shallow,1)...
                ./nansum(thickness_shallow,1);
        otherwise
            avg_shallow = nanmean(M_shallow,1);
    end
%     figure
%     for i = 1 ... :100 ... floor(length(avg_shallow)/3) ...
%             :length(avg_shallow)
%         subplot(1,2,1)
%         plot([0 20] , avg_shallow(i)*[1 1])
%         hold on 
%         plot(depth_alt(:,i),M_shallow(:,i))
%         hold off
%         title(sprintf('%i',i))
%         ylim([0 800])
%         subplot(1,2,2)
%         plot(1:i,avg_shallow(1:i))
%                 ylim([450 550])
%                 xlim([0 length(avg_shallow)])
%         pause(0.01)
%     end
end