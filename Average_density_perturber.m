%% Takes in 89DP, reshapes then perturbs density in a 10cm strip by an amount 'r', outputs new csv file


function  [] = Average_density_perturber(DP89_T20_densities,DP89_T20_depths)

r = 1.0:0.1:2;
%DP89_T20_new_densities= zeros();
Newcsv = reshape(DP89_T20_densities,100,20);
for i = 1:length(r)
    Altered = Newcsv(8:18,:)*r(i);
Newcsv(8:18,:) = Altered;
DP89_T20_new_densities = reshape(Newcsv,2000,1);
A(i) = [DP89_T20_depths,DP89_T20_new_densities];
writematrix(A,'Newcsvfile.csv')

end

%% plot check for perturbation

% subplot(1,2,1)
% plot(DP89_T20_densities,DP89_T20_depths)
% set(gca,'Ydir','reverse')
% xlabel('Density (kg m^{-3})')
% ylabel('Depth(m)')
% title('Original Density profile Top 20')
% hold on
% subplot(1,2,2)
% plot(A(:,2),A(:,1))
% set(gca,'Ydir','reverse')
% title('Perturbed profile')
% xlabel('Density (kg m^{-3})')
% ylabel('Depth(m)')
% legend('Original','Perturbed')
% saveas(gcf,'Original and Perturbed Density profiles','png')