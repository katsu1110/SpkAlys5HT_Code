function [minv] = nanmin(v)

v(v==0) = [];
if ~isempty(v)
    minv = min(v);
else
    minv = nan;
end