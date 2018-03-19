function [maxv] = nanmax(v)

v(v==0) = [];
if ~isempty(v)
    maxv = max(v);
else
    maxv = nan;
end
