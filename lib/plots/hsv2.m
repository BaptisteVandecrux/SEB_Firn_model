function map = hsv2(m)
%HSV    Hue-saturation-value color map
%   HSV(M) returns an M-by-3 matrix containing an HSV colormap.
%   HSV, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   An HSV colormap varies the hue component of the hue-saturation-value
%   color model.  The colors begin with red, pass through yellow, green,
%   cyan, blue, magenta, and return to red.  The map is particularly
%   useful for displaying periodic functions.  
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(hsv)
%
%   See also GRAY, HOT, COOL, BONE, COPPER, PINK, FLAG, PRISM, JET,
%   COLORMAP, RGBPLOT, HSV2RGB, RGB2HSV.

%   See Alvy Ray Smith, Color Gamut Transform Pairs, SIGGRAPH '78.
%   C. B. Moler, 8-17-86, 5-10-91, 8-19-92, 2-19-93.
%   Copyright 1984-2004 The MathWorks, Inc.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
m1 = m;
m=m+floor(m/2.7);
h = (0:m-1)'/max(m,1);
if isempty(h)
  map = [];
else
  map = hsv2rgb([h ones(m,2)]);
  map = map(end-m1+1:end,:);
%   map = flipud(map);
end
