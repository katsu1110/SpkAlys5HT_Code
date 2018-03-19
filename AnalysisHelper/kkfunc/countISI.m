function [isic] = countISI(exinfo)

len_exp = length(exinfo);
isic = nan(1, len_exp);
for i = 1:len_exp
    isic(i) = 100*sum(exinfo(i).isi > 0.02 & exinfo(i).isi < 2)/length(exinfo(i).isi);
end

close all
plot(1:len_exp, isic)
xlim([0.5 len_exp+0.5])
xlabel('experiments')
ylabel('ISI fraction (%): 0.02 < ISI < 2 ')
set(gca, 'box', 'off')
set(gca, 'TickDir', 'out')
axis square