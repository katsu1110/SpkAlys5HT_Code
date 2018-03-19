function m = nonzeromedian(v)

m = median(v(~ismember(v, 0) & ~isnan(v)));