function [depth_hor] = track_horizon(depth_start,ind_start,step,...
                            depth_act, compaction, H_surf,ind_stop)


length_out = size(depth_act,2);
if ~exist('ind_stop','var')
    ind_stop = length_out;
end
depth_hor = NaN(length_out,1);
depth_hor(ind_start) = depth_start;

for i = (ind_start+step):step:ind_stop
    depth_hor(i) = max(0,depth_hor(i-step) + ...
       (H_surf(i)-H_surf(i-step)));
%     depth_hor(i) = depth_hor(i-step) ;
    depth_mod = depth_act(:,i);
    comp_mod = compaction(:,i)*step;
    depth_new = [depth_mod; depth_hor(i)];
    [depth_new, ind_sorted] = sort(depth_new);
    % index of the layer containing horizon in depth_new
    ind_hor = find(ind_sorted == length(depth_mod)+1);
    [comp_new] = ConvertToGivenDepthScale(depth_mod, comp_mod, ...
        depth_new,'extensive');

    comp_tot = sum(comp_new(ind_hor+1:end));
    depth_hor(i) = depth_hor(i) + comp_tot;
end

% interpolating between the steps
 depth_hor(isnan(depth_hor))=interp1(find(~isnan(depth_hor)),...
    depth_hor(~isnan(depth_hor)),...
    find(isnan(depth_hor)),'linear'); 

end
