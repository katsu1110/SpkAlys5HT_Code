
function mu = getMean(ax)
% get mean from title

mu(1) = getMedianHelper(ax.Title.String{2});
mu(2) = getMedianHelper(ax.Title.String{3});
end

function mu = getMeanHelper(s)
[~, idx] = regexp(s, '\mu=');

if strcmp(s_5HT(idx+1), '-')
    mu = str2double( s_5HT(idx+1:idx+6) );
else
    mu = str2double( s_5HT(idx+1:idx+5) );
end
end