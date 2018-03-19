function [r, p] = nancorr(x,y, name)

nanidx = union(find(isnan(x)), find(isnan(y)));
x(nanidx) = [];
y(nanidx) = [];

try
    [r,p] = corr(x,y,'type',name);
catch
    r = nan; p = nan;
end