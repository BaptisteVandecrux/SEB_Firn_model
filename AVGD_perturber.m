%% Takes in 89DP, reshapes then perturbs density in a 10cm strip by an amount 'r', outputs new csv file
% Average_density_perturber

DP_1989 = readtable('DensityProfile_KAN-U_1989_metres.csv');
Original_DP89_T20_depths = table2array(DP_1989(:,1));
Original_DP89_T20_densities = table2array(DP_1989(:,2));

DP89_T20_depths = table2array(readtable('DP89_Top20_depths'));
DP89_T20_densities = table2array(readtable('DP89_Top20_densities'));
rr = table2array(readtable('r.txt'));
r = rr*2/100 + 1;
Newcsv = reshape(DP89_T20_densities,100,20);
Altered = Newcsv(8:18,:)*rr;
Newcsv(8:18,:) = Altered;
DP89_T20_new_densities = reshape(Newcsv,2000,1);

Original_DP89_T20_densities(1:2000) = DP89_T20_densities;
Perturbed_Densities = Original_DP89_T20_densities;
A = [Original_DP89_T20_depths,Perturbed_Densities];
writematrix(A,'./Input/Initial state/density/Site_J_perturbations.csv')
% 
% %subplot(1,2,1)
% plot(Original_DP89_T20_densities,Original_DP89_T20_depths)
% set(gca,'Ydir','reverse')
% xlabel('Density (kg m^{-3})')
% ylabel('Depth(m)')
% %ylim([0 200])
% title('Original Density profile Top 20')
% hold on
% %subplot(1,2,2)
% plot(A(:,2),A(:,1))
% set(gca,'Ydir','reverse')
% title('Perturbed profile')
% xlabel('Density (kg m^{-3})')
% ylabel('Depth(m)')
% %ylim([0 200])
% %legend('Original','Perturbed')
% saveas(gcf,'Original and Perturbed Density profiles','png')