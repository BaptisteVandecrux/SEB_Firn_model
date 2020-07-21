function [data_out] = ChangeOrigin(data_in)
% in input:
% 0: main station
% 1: CP2 at 6.5 km
% 2: Swiss Camp at 97 km
% 3: KAN_U at 67 km
% 4: RCM
% 5: MODIS
% 9: NOAA tower at 2 km
% 10: Miller's station at 2 km
%
% in output:
% 0: main station
% 0<X<100: adjusted from secondary AWS located X km away from GC-Net station  
% 101: adjusted from RACMO2.3p2
% 102: calculated from SRin and nearest daily MODIS albedo (or MODIS albedo
% climatology if before 2000)

data_out = data_in;
data_out(data_in ==4) = 101;
data_out(data_in ==5) = 102;
data_out(data_in ==1) = 6.5;
data_out(data_in ==2) = 97;
data_out(data_in ==3) = 67;
data_out(data_in ==9) = 2;
data_out(data_in ==10) = 2;
end