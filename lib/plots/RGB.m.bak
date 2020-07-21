function RGB= RGB(index_or_name)
%
% OVERVIEW
%
% The Matlab function RGB() converts a color index (whole number from
% 1-21), English name of a color (string), or RGB triple with whole number
% components in {0, 1, ..., 255} into an RGB triple with real-valued
% components in [0,1].  RGB() allows the user to access a set of 21 colors
% via their common English names.  For eight of these colors that have more
% than one common name, the program accepts alternative names, reducing the
% memory burden on the user.
%
% INPUTS
%
% - index_or_name can be a color index (whole number from 1-21), a string
% containing the name of a color in lower case, an RGB triple with elements
% in [0,1], or an RGB triple with elements in [0,255].  Note that some
% colors have more than one name, in which case any of these may be used.
% See the code for the list of color names.
%
% OUTPUTS
%
% - RGB is a length-3 vector of RGB components that can be used as a color
% specifier with any of the Matlab graphics functions.
%
% If the input is an RGB triple with elements in [0,1], it is returned to
% the calling program without modification.
%
% If the input is an RGB triple with elements in [0,255], it is scaled by
% 1/255 and then returned.
%
% If the input is a color index (1-21), it is converted to an RGB triple
% via direct table lookup.
%
% If the input is the name of a color, a search is done to find a matching
% name, and the corresponding RGB triple is returned.
%
%
% DEPENDENCIES
%
% RGB.m depends on cell2mat.m and lookup.m.
%
%
% FUTURE ENHANCEMENTS
%
% I plan to at some point add code to allow the user to define additional
% colors and change the order of the colors.
%
% Dr. Phillip M. Feldman 9 April 2009


% Section 1: Define names and RGB values of 21 colors.  Note that some
% colors have more than one name.

colors= {
1    0    0    'red';
194/255 78/255 67/255 'rouge';
194/200 78/200 67/200 'rouge clair';
1    0.6  0    'orange';
0.25 0.25 0.9  {'blue','cobalt blue'};
0    1    0    {'green','light green'};
63/255 182/255 83/255    {'vert'};
0.2    0    0    'dark red';
1    0.5    0.5    'light red';
1    0.8    0.8    'light light red';
0    1    1    {'cyan','baby blue'};
0.5  0.5  1    'light blue';
0.8  0.8  1    'light light blue';
0.8  0.5  0    'brown';
0.9  0.9  0.9  'light light gray';
0.7  0.7  0.7  'light gray';
0.6  0.6  0.6  'gray';
0.5  0.5  0.5  'dark gray';
1    1    0.6  'cream';
0    0.5  0    {'dark green','forest green'};
0    0.25  0    {'dark dark green'};
0.8  0.8  0.8  {'gray','light gray'};
0    0    0    'black';
0.5    1    0.5    {'light light green'};
1    0.5  0.5  'peach';
1    1    0    'yellow';
223/255 199/255 79/255    'jaune';
0    0    0.8  {'dark blue','navy blue'};
0.1    0.1    0.3  {'dark dark blue'};
0.8  0.8  0.8  {'gray','light gray'};
0.5  0    0.9  'purple';
1    0.5  1    {'magenta','pink'};
1    0.8  1    {'light pink'};
0    0.8  0.8  {'aqua','turquoise'};
0.9  0.75 0    'gold';
1    1    1    'white';
0.9             0.7         1 'light purple';
};

% Extract names from rightmost column of colors cell array:
names= colors(:,4);

% Strip off rightmost column of colors cell array and convert remaining
% columns to a matrix:
RGBs= cell2mat(colors(:,1:3));


% Section 2: Convert index_or_name to an RGB triple.

if nargin ~= 1
   error('This function must be called with exactly one argument.');
end

if isnumeric(index_or_name)
   index= index_or_name;


   if length(index)==3 & all(index>=0)

      % If contents of index_or_name are an RGB triple with elements in
      % [0,1], return them as output without modification:

      if all(index<=1)
         RGB= index;
         return
      end

      % If contents of index_or_name are an RGB triple with elements in
      % [0,255], scale by 1/255 and return as output:

      if all(index<=255)
         RGB= index/255;
         return
      end

   end

   if length(index) > 1
      error('When calling with a color index, specify a single number.');
   end
   if ismember(index,1:21)
      RGB= RGBs(index,:);
      return
   end
   error('A color index must be a whole number between 1 and 21.');
end

if isa(index_or_name,'char')
   index = lookup(names,index_or_name);
   if index
      RGB = RGBs(index,:);
   else
      fprintf(2, ['Warning: Unknown color name "%s".  ' ...
        'Substituting black.\n'], index_or_name);
      RGB= [0 0 0];
   end
   return
end

error('Input argument has unexpected data type.');
