function n= lookup(c, str)
%
% OVERVIEW
%
% lookup(c, str) compares the string in str against each element of the
% cell array c in turn.  If str is found in any element of c, the index of
% the first such element is returned.  Otherwise, the value zero is
% returned.
%
% INPUTS
%
% c: A one-dimensional cell-array.  Each element of this cell array is
% either a character string, or another one-dimensional cell-array
% containing character strings.  Note: Deeper nesting of cell arrays is not
% supported.
%
% EXAMPLES
%
% All three of the following function calls return the value 2:

% lookup({'cat','dog','fish'},'dog')
% lookup({{'c','cat'},{'d','dog'},'fish'},'dog')
% lookup({{'c','cat'},{'d','dog'},'fish'},'d')
%
% Dr. Phillip M. Feldman
% 3 June 2008

if (nargin ~= 2)
   error('This function requires exactly two calling arguments.');
end
if (~iscell(c))
   error('The first calling argument must be a cell array.');
end
if (~isa(str,'char'))
   error('The second calling argument must be a character string.');
end

for i= 1 : length(c)
   ndx(i,1)= any(strcmpi(c{i},str));
end

n= find(ndx, 1, 'first');

% We can't return an empty matrix, so convert an empty matrix to a zero:

if (isempty(n)), n= 0; end

end
