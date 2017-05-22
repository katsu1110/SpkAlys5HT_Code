function med = getMedian(ax)
% get median from title

med(1) = getMedianHelper(ax.Title.String{2});
med(2) = getMedianHelper(ax.Title.String{3});
end

function med = getMedianHelper(s)
[~, idx] = regexp(s, 'med=');

if strcmp(s(idx+1), '-')
    med = str2double( s(idx+1:idx+6) );
else
    med = str2double( s(idx+1:idx+5) );
end
end