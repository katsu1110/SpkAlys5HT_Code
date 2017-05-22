
function stars = getSigStars(p)

if p < 0.001
    stars = '***';
elseif p<0.01
    stars = '**';
elseif p<0.05
    stars = '*';
elseif p<0.08
    stars = 'o';
else
    stars = '';
end
    
end
