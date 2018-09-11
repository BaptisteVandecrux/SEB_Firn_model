function [smoothedarray] = MySmooth(array, span)
% ind_bin = discretize(array,[0 common_depth_scale]);
% smoothedarray = accumarray(ind_bin,array,[],@nanmean);

smoothedarray=array;
% old version 1
% for i=1:length(array)
%     ind = (i-span):(i+span);
%     ind = ind(ind > 0);
%     ind = ind(ind < length(array));
%     smoothedarray(i) = nanmean(array(ind));
% end
% old version 2
for i=1:floor(length(array)/span)+1
    ind = ((i*span+1):((i+1)*span))-span;
    ind = ind(ind > 0);
    ind = ind(ind <= length(array));
    smoothedarray(ind) = nanmean(array(ind));
end
smoothedarray(isnan(array))=NaN;
end
