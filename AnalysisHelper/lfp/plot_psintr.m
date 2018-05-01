function plot_psintr(psintr)
% plot results from 'pupil_interaction.m'
% for RC analysis: tuning curve, type II regression
% DRUG, PS, sps: base vs drug, lps: base vs drug

close all;
h = figure;
switch psintr.drugname
    case '5HT'
        dcol = 'r';
    otherwise
        dcol = 'k';
end
name = {{'base', 'drug','drug'}, {'lps','sps','ps'}, {3,1,'sps_drug'},{4,2,'lps_drug'}};
yall = [];
for i = 1:4    
    % tuning curve
    if i < 3
        xv1 = [psintr.(['rcsub_' name{i}{1}]).results.stm.val];
        yv1 = sqrt([psintr.(['rcsub_' name{i}{1}]).results.stm.peak]);
        xv2 = [psintr.(['rcsub_' name{i}{2}]).results.stm.val];
        yv2 = sqrt([psintr.(['rcsub_' name{i}{2}]).results.stm.peak]);
    else
        xv1 = [psintr.rcsub(name{i}{1}).results.stm.val];
        yv1 = sqrt([psintr.rcsub(name{i}{1}).results.stm.peak]);
        xv2 = [psintr.rcsub(name{i}{2}).results.stm.val];
        yv2 = sqrt([psintr.rcsub(name{i}{2}).results.stm.peak]);
    end
    yall = [yall, yv1, yv2];
    subplot(2,4,i)
    l = zeros(1,2);
    l(1) = plot(xv1, yv1, '-ok');
    hold on;
    l(2) = plot(xv2, yv2, '-or');
    xlabel('orientation (^o)')
    ylabel('mean firing rate (spk/s)')
    if i==2
        legend(l, name{i}{1}, name{i}{2}, 'location', 'northwest')
    else
        legend(l, name{1}{1}, name{1}{2}, 'location', 'northwest')
    end
    legend('boxoff')
    title(name{i}{3})
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    % type II regression
    subplot(2,4,i+4)
    m = yv1;
    m = max(m);
    xv = yv1/m;
    yv = yv2/m;
    scatter(xv, yv, 50, 's', ...
        'markerfacecolor', dcol, 'markerfacealpha', 0.4, ...
        'markeredgecolor', 'w', 'markeredgealpha', 0.8)
    hold on;
    minima = min([xv, yv])*0.95;
    maxima = max([xv, yv])*1.05;
    plot([minima, maxima],[minima, maxima],'--k')
    hold on;
    reg = psintr.type2reg.(name{i}{3});
    plot([minima, maxima], reg(2)*[minima, maxima]+reg(1), ...
        '-', 'color', dcol, 'linewidth',1)
    if i==2
        xlabel(name{i}{1})
        ylabel(name{i}{2})
    else
        xlabel(name{1}{1})
        ylabel(name{1}{2})
    end    
    axis([minima maxima minima maxima])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end
yrange = [min(yall)*0.95, max(yall)*1.05];
for i = 1:4
    subplot(2,4,i)
    ylim(yrange)
end