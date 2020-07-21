function [h] = stairs0(x,y, varargin)
    if isrow(x)
        x = x';
    end
    if isrow(y)
        y = y';
    end
    x = [x; x(end) + (x(end)-x(end-1))];
    y = [y; y(end)];
    h = stairs(x,y, varargin{:});
end