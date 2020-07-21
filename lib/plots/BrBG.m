function h = BrBG(m)
% blues colormap from brewermap

if nargin < 1, m = size(get(gcf,'colormap'),1); end
[h,~ ,~] = brewermap(m,'BrBG');
h = flipud(h);
