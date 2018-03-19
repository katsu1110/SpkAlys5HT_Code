function m = nonzerostd(v)

m = std(v(~ismember(v, 0) & ~isnan(v)));