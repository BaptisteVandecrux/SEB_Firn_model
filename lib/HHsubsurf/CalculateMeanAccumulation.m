function c = CalculateMeanAccumulation(time,snowfall_weq, c)

% here we calculate the accumulation rate (m_weq/year) that is used in dry densfication
years = unique(floor(time));
precip_year = NaN(length(years)-2,1);
% only calculating total precipitation on full years
for i = 2:length(years)-1
    precip_year(i-1) = sum(snowfall_weq(floor(time)== years(i)));
end
c.mean_accum = mean(precip_year);

% if the period does not contain one full year, the accumulation rate is
% then calculated from mean accumulation rate (m_weq/hr) 
if isnan(c.mean_accum)
    c.mean_accum = sum(snowfall_weq)/length(snowfall_weq)*24*365;
end

end