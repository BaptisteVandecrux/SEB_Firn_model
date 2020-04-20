function [varargout] = InterpolateBackToOriginal(M, varargin)
        ipfac2 = 1:(M+1);
        for i = 1:length(varargin)
            varargin{i}        = interp1(varargin{i},ipfac2);
        end
        varargout = varargin;
end