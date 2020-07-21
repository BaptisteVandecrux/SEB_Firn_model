function [f1, f2, EventTime, EventDepth, EventMagnitude] = PipingStudy(time_obs,T_obs, depth_thermistor, ...
    Surface_Height, varargin)

EventTime=[];
EventDepth=[];
EventMagnitude=[];

% default parameters
param.dt = 1; %time step in hour
param.Threshold = 1; %threshold for the peak detection
param.Range = -25:0; %range for the plot
param.Visible = 'on';
param.Title = '';
%assigning user defined parameters

if floor(length(varargin)/2) ~= (length(varargin))/2
    fprintf('Error: Wrong number of arguments\n')
    return
else
for i = 1:2:length(varargin)
    param.(varargin{i}) = varargin{i+1};
end
end

T_obs(depth_thermistor==0)=NaN;
depth_thermistor(depth_thermistor==0)=NaN;
T_obs(T_obs==-999)=NaN;

TT_obs = repmat(time_obs',size(T_obs,2),1);
% percol = imregionalmax(T_obs);
i_percol = [];
j_percol = [];
for j = 2:size(T_obs,2)-1
    [~, lcs] = findpeaks(T_obs(:,j),'MinPeakProminence',param.Threshold);
    for i=1:length(lcs)
        if ~isnan(T_obs(lcs(i),j-1))&&~isnan(T_obs(lcs(i),j+1))
            if T_obs(lcs(i),j)>T_obs(lcs(i),j+1) &&...
                    T_obs(lcs(i),j) > T_obs(lcs(i),j-1)
                if T_obs(lcs(i),j-1)<-0.2
                    j_percol = [j_percol; j];
                    i_percol = [i_percol; lcs(i)];
                end
            end
        end
    end
end

T_obs(T_obs==-999) = NaN;
T_obs(depth_thermistor<0)=NaN;
depth_thermistor(depth_thermistor<=0)=NaN;
depth_thermistor(isnan(T_obs))=NaN;

for i = 1:size(T_obs,2)
    depth_thermistor(:,i) = depth_thermistor(:,i) - Surface_Height;
end

i_percol_old = i_percol;
j_percol_old = j_percol;

ind_bad = [];
for i = 1:length(i_percol)
    if sum(isnan(T_obs(i_percol(i),1:(j_percol(i)-1))))==length(1:(j_percol(i)-1))
        ind_bad = [ind_bad, i];
    end
end
i_percol(ind_bad) = [];
j_percol(ind_bad) = [];

% Removing points that have nan temp before or after
ind_bad = [];
for j = 1:length(j_percol)
    if isnan(T_obs(i_percol(j)-1,j_percol(j))) || ...
            isnan(T_obs(i_percol(j)+1,j_percol(j)))
        ind_bad = [ind_bad, j];
    end
end
i_percol(ind_bad) = [];
j_percol(ind_bad) = [];

f1 = figure('Visible',param.Visible);
PlotTemp(TT_obs, depth_thermistor', T_obs',...
    'PlotTherm', 'yes',...
    'PlotIsoTherm', 'yes',...
    'Range', param.Range);
ylim([-1, 10])
title(param.Title)

plot(time_obs,-Surface_Height,'LineWidth',2)

% scatter(time_obs(i_percol_old), ...
%     depth_thermistor(sub2ind(size(depth_thermistor),i_percol_old,j_percol_old)), 50,...
%     'filled', 'r')

scatter(time_obs(i_percol), ...
    depth_thermistor(sub2ind(size(depth_thermistor),i_percol,j_percol)), 50,...
    'filled','k')
text(time_obs(i_percol)'+0.5, ...
    depth_thermistor(sub2ind(size(depth_thermistor),i_percol,j_percol))'+0.5,...
    cellstr(num2str([1:length(i_percol)]')))

%% recalculate periods based on new percolation events
ind_start=[];
ind_end = [];
% clear ind_start ind_end

for i = 1:length(i_percol)
    ii = i_percol(i);
    while T_obs(ii,j_percol(i))>T_obs(ii,j_percol(i)-1) && ii > 1
    ii = ii-1;
    end
    ind_start(i) = ii;
    ii = i_percol(i);
    while T_obs(ii,j_percol(i))>T_obs(ii,j_percol(i)-1)&& ii < size(T_obs,1)
    ii = ii+1;
    end
    ind_end(i) = ii;
end
percol_event = unique([ind_start' ind_end' j_percol],'rows');

%% Plotting events
i_new = 1;  

f2{1} = figure('Visible',param.Visible);
if ~isempty(i_percol)
[ha, pos] = tight_subplot(3,4, [0.07, 0.025], [0.06, 0.04],[0.02, 0.02]);

for i = 1: length(i_percol) %size(percol_event,1)
        
    time_ind_prev = max(i_percol(i) - 3*24/param.dt, 1);
    time_ind_next = min(i_percol(i) + 3*24/param.dt, length(time_obs));
    EventTime = [EventTime, time_obs(i_percol(i))];
    EventDepth = [EventDepth, depth_thermistor(i_percol(i),j_percol(i))];
    EventMagnitude = [EventMagnitude, T_obs(i_percol(i),j_percol(i)) - ...
        (T_obs(time_ind_prev,j_percol(i)) +T_obs(time_ind_next,j_percol(i)))/2];

    if i_new==13
        legend('event', '-1 day', '+1 day','Location','SouthWest')
        set(gcf,'Visible',param.Visible)
        f2{length(f2)+1} = figure('Visible',param.Visible);
        i_new=1;
        [ha, pos] = tight_subplot(3,4, [0.07, 0.025], [0.06, 0.04],[0.02, 0.02]);
    end

    axes(ha(i_new))
    plot(T_obs(i_percol(i), :),depth_thermistor(i_percol(i), :),'k', 'LineWidth',2)
    hold on
    
    plot(T_obs(time_ind_prev, :),depth_thermistor(time_ind_prev, :),'r')
    plot(T_obs(time_ind_next, :),depth_thermistor(time_ind_next, :),'b')

    scatter(T_obs(i_percol(i),j_percol(i)), depth_thermistor(i_percol(i),j_percol(i)),70,'filled')
    set(gca,'Ydir','reverse')
    axis tight
    [leg,obj,~,~] = legend('','Location','SouthWest');
    p = leg.Position;
        title(datestr(time_obs(i_percol(i))))

    legend('off')
    %     BoxPosition=[p(1)+0.133 p(2)+0.025 0  p(4)];
    annotation('textbox',...
        'Units','normalized',...
        'Position',p,...
        'String',num2str(i),...
        'FitBoxToText','on',...
        'EdgeColor','k',...
        'FontWeight','bold',...
        'FontSize',10.5,...
    'Interpreter','none');
    xlim([-15 0]);
    ylim([0 10]);
    i_new=i_new+1;

end
legend('event', '-1 day', '+1 day','Location','SouthEast')

while i_new < 13
    axes(ha(i_new))   
    set(gca, 'Visible','off')
    i_new = i_new+1;
end
end
