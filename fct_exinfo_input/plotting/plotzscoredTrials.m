function plotzscoredTrials(exinfo)
% plot the time course of spiking activity by using the zscored spike
% counts in the first 70 trials.



idx =   find( cellfun(@isempty, strfind({exinfo.fname}, 'all.grating')) & ...
                cellfun(@isempty, strfind({exinfo.fname_drug}, 'all.grating')) );
exinfo = exinfo(idx);

exinfo = exinfo(~[exinfo.isRC]);

ntrials = 70;
t = 1:ntrials;

figure;
subplot(2,1,1);

for kk = 1:length(exinfo)

    temp = exinfo(kk).trials_c1.zspkCount;
    temp(length(temp):300) = nan;
    temp = temp(1:ntrials);
    z_base(kk,:) = temp;

    plot(t, temp, 'k', 'LineWidth', 0.5);hold on;
end


mn = nanmean(z_base);
sd = nanstd(z_base);


fill([t fliplr(t)], [mn+sd fliplr(mn-sd)], 'r', 'FaceAlpha', .3, 'EdgeColor', 'none'); hold on;
plot(t, mn, 'r', 'LineWidth', 2);
title('baseline');
ylim([-5 5])



%%%
subplot(2,1,2);

exinfo = exinfo([exinfo.is5HT] & [exinfo.nonparam_ratio]<1);
for kk = 1:length(exinfo)
    temp = exinfo(kk).trials_c1.zspkCount;
    temp(length(temp):300) = nan;
    temp = temp(1:ntrials);
    z_drug(kk,:) = temp;
    
    plot(t, temp, 'k', 'LineWidth', 0.5); hold on;

end

mn = nanmean(z_drug);
sd = nanstd(z_drug);
t = 1:length(mn);

fill([t fliplr(t)], [mn+sd fliplr(mn-sd)], 'r', 'FaceAlpha', .3, 'EdgeColor', 'none'); hold on;
plot(t, mn, 'r', 'LineWidth', 2);
title('5HT')
ylim([-5 5])
return 
end

