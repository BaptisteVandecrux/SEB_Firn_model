function h = blues(m)
% blues colormap from brewermap

if nargin < 1, m = size(get(gcf,'colormap'),1); end
[h,~ ,~] = brewermap(m,'blues');
h = flipud(h);
