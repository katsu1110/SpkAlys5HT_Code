function m = nonzeromean(v)

m = mean(v(~ismember(v, 0) & ~isnan(v)));