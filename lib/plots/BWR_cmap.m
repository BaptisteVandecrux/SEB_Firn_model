function cmap = BWR_cmap(num_color)
    if nargin < 1
       m = size(get(gcf,'colormap'),1);
    end
    % defining the colormap
    start_color = [0,0, 1];
    end_color = [1, 0, 0];

    cmap = ones(num_color,3);
    for i = 1:floor(num_color/2)
        cmap(i,:) = start_color + ...
            ([1, 1, 1] - start_color) * (i-1)/floor(num_color/2);
        cmap(num_color-i+1,:) = end_color ...
            +([1, 1, 1] - end_color) * (i-1)/floor(num_color/2);
    end
        cmap(i,:) = [1 1 1];
        cmap(num_color-i+1,:) = [1 1 1];
    cmap(end,:) = [0.8 0 0];
end
