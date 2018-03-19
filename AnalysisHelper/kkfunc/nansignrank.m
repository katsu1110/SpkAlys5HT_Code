function [p,h,stats] = nansignrank(x,y)

try
    [p,h,stats] = signrank(x,y);
catch
    p = nan; h = nan; stats.zval = nan; stats.signedrank = nan;
end