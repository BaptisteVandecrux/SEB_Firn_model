function [x,y] = temporal_shift(x_in,y_in,varname, GapFillMethod)
% t = data_AWS_GITS.time(ind_GITS);
% x = data_AWS_GITS.(VarList2{k})(ind_GITS);
% y=data_surf{ii}.(VarList{k})(ind_mod_GITS)*fac;
    x = x_in;
    y = y_in;

    % double check regarding thhe synchronization of the two time series
        x_scaled = (x - nanmean(x)) / nanstd(x);
%         x_scaled(isnan(x)) = 0;
        y_scaled = (y - nanmean(y)) / nanstd(y);
%         y_scaled(isnan(y)) = 0;
        isOK=isfinite(x_scaled) & isfinite(y_scaled);   % both rows finite (neither NaN)
        [acor, lag] = xcorr(x_scaled(isOK),y_scaled(isOK),20);
        
%         figure
%         subplot(2,1,1)
%         hold on
%         plot(y_scaled)
%         plot(x_scaled)
% 
%         subplot(2,1,2)
%         plot(lag,acor)
        if find(lag(acor==max(acor))) ~= 0
            shift = lag(acor==max(acor));
            disp(['Lag detected: ', varname ' ' num2str(shift) ' hr'])

            if GapFillMethod == 1 
                if shift >0
                    y(shift+1:end) = y(1:end-shift);
                else
                    y(1:end+shift) = y(-shift+1:end);
                end
            else 
                shift = -shift;
                if shift >0
                    x(shift+1:end) = x(1:end-shift);
                else
                    x(1:end+shift) = x(-shift+1:end);
                end
%                 x_scaled = (x - nanmean(x)) / nanstd(x);
%                 y_scaled = (y - nanmean(y)) / nanstd(y);
%                 isOK=isfinite(x_scaled) & isfinite(y_scaled);   % both rows finite (neither NaN)
%                 [acor, lag] = xcorr(x_scaled(isOK),y_scaled(isOK),20);
            end
%             subplot(2,1,1)
%             hold on
%             plot(x_scaled)
%             subplot(2,1,2)
%             hold on
%             plot(lag,acor)
        end
end