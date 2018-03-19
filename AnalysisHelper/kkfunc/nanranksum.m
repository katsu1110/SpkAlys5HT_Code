function [p,h,stats] = nanranksum(x,y)

try
    [p,h,stats] = ranksum(x,y);
catch
    p = nan; h = nan; stats.zval = nan; stats.ranksum = nan;
end