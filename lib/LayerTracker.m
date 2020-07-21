

% depth_act_2 = depth_act;
% midpoint_depth = depth_act(2:end,:) - thickness_act./2;
% midpoint_depth_act_2 = vertcat(zeros(size(depth_act_2(1,:))), midpoint_depth);
% 
rho_2 = vertcat(rho_all, rho_all(end,:));
% rho_3 = vertcat(rho, rho(end,:));
% slwc_2 = vertcat(slwc, slwc(end,:));
% snic_2 = vertcat(snic(1,:),snic);
TT_2 = vertcat(TT,TT(end,:));
% for i = 1:size(depth_act_2,1)
%     depth_act_2(i,:) = depth_act_2(i,:) - H_surf' + H_surf(1);
%     midpoint_depth_act_2(i,:) = midpoint_depth_act_2(i,:) - H_surf' + H_surf(1);
% end
tic
f = figure('Visible', vis);
PlotTemp(TT_2,depth_act,rho_2,...
    'Range',300:50:900,...
    'PlotTherm','no',...
    'PlotIsoTherm','no',...
    'Interp','off');
hold on
[C, h] = contour(TT_2(:, 1:end), depth_act(:, 1:end), ...
    repmat([0:size(TT_2(:, 1:end),1)-1]',1,size(TT_2(:, 1:end),2)),0:201,'LineWidth',0.001,'Color','k');
% clabel(C,h);
axis tight
toc
set(f,'Visible','on')

% xlabel('Time')
% ylabel('Real depth (m)') 
ylim([min(-1,min(-H_surf)), 10])
% ylim([50, 80])
print(f,sprintf('%s/LayerTracker',c.OutputFolder),'-dtiff')