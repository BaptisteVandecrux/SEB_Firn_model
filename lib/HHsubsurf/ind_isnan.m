function [i] = ind_isnan(T,str)
if strcmp(str, 'beginning')
    i=1;
    while isnan(T(i))
        i=i+1;
    end
elseif strcmp(str, 'end')
    i=length(T);
    while isnan(T(i))
        i=i-1;
    end
else
    dips('ERROR: string not valid')
end
    
end