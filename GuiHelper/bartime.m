function bartime( dat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure;

%%%
binranges = exp(-2.1:0.1:1.6); 

for i = 1:max(dat.x)-2
    [nelements(i, :), center] = hist(dat.y(dat.x==i), binranges);
    nelements(i, :) = nelements(i, :) / max(dat.y(dat.x==i));
    posy(i) = find( median(dat.y(dat.x==i)) < center, 1, 'first');
end

%%%

d = conv2(nelements', ones(3,1));

for j=1:max(dat.x)-2
   d(:, j) = d(:, j) ./ max(d(:, j));   
end

imagesc(d); hold on; plot(1:6, posy, 'w', 'LineWidth', 2');
set(gca, 'YDir', 'normal', 'YTick', 1:2:length(binranges), 'YTickLabel', cellstr(num2str(binranges(1:2:end)')) );


end