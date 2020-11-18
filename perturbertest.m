Newcsv = reshape(DP89_T20_densities,100,20);
Altered = Newcsv(8:18,:)*2;
Newcsv(8:18,:) = Altered;
DP89_T20_new_densities = reshape(Newcsv,2000,1);
A = [DP89_T20_depths,DP89_T20_new_densities];
writematrix(A,'P10.csv')
