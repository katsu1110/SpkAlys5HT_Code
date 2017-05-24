function exinfo = phasePlot( exinfo, ex0, ex2 )
% stimulus phase as time vs psth
% this plot reflects evidence for simple or complex behavior


hnew = figure('Name', exinfo.figname, 'UserData', exinfo, ...
    'Position', [288   144   849   853]);

%%
% load baseline and 5HT/NaCl plots
getPhaseSelectivity(ex0, 'stim', exinfo.param1, 'plot',true);
h0 = gcf;
ax0 = h0.Children(2:end);

getPhaseSelectivity(ex2, 'stim', exinfo.param1, 'plot',true);
h2 = gcf;
ax2 = h2.Children(2:end);



%% Plot in new Figure
figure(hnew)
% psth with stimulus
for i = 1:length(ax0)-2
    s0(i) = subplot(length(ax0)-2, 3, (i*3)-2);
    copyobj(ax0(i+2).Children, s0(i));
    title(ax0(i+2).Title.String)
    x = ax0(i+2).Children.XData;
    xlim([x(1) x(end)]); xlabel spk/s;
end
for i = 1:length(ax2)-2
    s2(i) = subplot(length(ax2)-2, 3, (i*3)-1);
    copyobj(ax2(i+2).Children, s2(i));
    title(ax2(i+2).Title.String); 
    xlim([x(1) x(end)]); xlabel spk/s;
end

s0(1).Title.String = sprintf(['Baseline \n' s0(1).Title.String]);
s2(1).Title.String = sprintf([exinfo.drugname '\n' s2(1).Title.String]);

% f1/f0 vs stimulus
s(1) = subplot(2, 3, 3);
plot(vertcat(ax0(2).Children.XData)', ...
    vertcat(ax0(2).Children.YData)'); hold on
plot(vertcat(ax2(2).Children.XData)', ...
    vertcat(ax2(2).Children.YData)', getCol(exinfo)); 
legend('baseline', exinfo.drugname);
xlabel stim;    ylabel f1/f0;

% f1/f0 vs stimulus - Baseline
s(2) = subplot(4, 3, 9);
plot(ax0(1).Children(1).XData, ax0(1).Children(1).YData, '--', 'Color', lines(1)); hold on;
plot(ax0(1).Children(1).XData, ax0(1).Children(2).YData, '-', 'Color', lines(1)); 
title(['Baseline: ' ax0(2).Title.String])
legend( 'f0', 'f1')


% f1/f0 vs stimulus - Drug
s(3) = subplot(4, 3, 12);
plot(ax2(1).Children(1).XData, ax2(1).Children(1).YData, '--', 'Color', getCol(exinfo)); hold on;
plot(ax2(1).Children(1).XData, ax2(1).Children(2).YData, '-', 'Color', getCol(exinfo)); 
title([exinfo.drugname ':' ax2(2).Title.String])
legend( 'f0', 'f1')

set(s, 'XTick', ax0(2).XTick, 'XTickLabel', ax0(2).XTickLabel, ...
    'XLim', ax0(2).XLim); 
set(findobj(hnew, 'Type', 'axes'), 'FontSize', 8)

xlabel stim;


%%
savefig(hnew, exinfo.fig_phase);

close(hnew); close(h0); close(h2);

end

