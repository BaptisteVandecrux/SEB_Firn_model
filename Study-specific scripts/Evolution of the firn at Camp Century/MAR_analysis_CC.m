    path = 'C:/Users/bav/OneDrive - Geological survey of Denmark and Greenland/Data/RCM/MAR/MARv3.9_rcp8.5/ACCESS1.3';

% finfo = ncinfo([path, '/' file]);
% CEN 77.1826N -61.1127E 1886 m a.s.l.
% # GITS polar stereo
% -386600 -1344500

years = 2006:2100;
smb = years;
runoff = smb;

for year = years
    year
    file = ['/MARv3.9-yearly-ACCESS1.3-rcp85-', num2str(year), '.nc'];
    RU = ncread([path, '/' file], 'RU');
    SMB = ncread([path, '/' file], 'SMB');
    time = ncread([path, '/' file], 'time');
    x = ncread([path, '/' file], 'x');
    y = ncread([path, '/' file], 'y');
    [~, i] = min (abs(x+386600));
    [~, j] = min (abs(y+1344500));    
    
    runoff(year-years(1)+1) = RU(i,j);
    smb(year-years(1)+1) = SMB(i,j);
    
%     figure
%     surf(SMB')
%     hold on
%    contour3(SMB','LevelList',[0],'LineColor' ,'k','LineWidth',2)
%      scatter3(i,j,SMB(i,j),'o','MarkerFaceColor','r')
%     hold 
%     shading interp
%     title(year)
%     colorbar
end

%%
figure
hold on
% scatter(years,runoff,100,'o','fill')
stairs(years,smb)
legend('SMB at Camp Century')
data_access = table(years',smb',runoff');
writetable(data_access,'MAR_ACCESS1.3_CC_smb_runoff.txt')

%% 
