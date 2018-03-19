function printfigureFFdiff
% this function loads the raw figures, processes the appearance
% and arranges them all in a new figure. the new figure is saved as pdf, if
% wanted

close all

[figPars, axPars] = setPlotPars;
% figPos = [10 10 21 29.7]; % this is in cm
figPos = [5 5 21 20]; % this is in cm
h = figure(figPars);
set(gcf, 'position', figPos, 'paperposition', figPos);

xbegin = 2.5;
ybegin = 15;
sq = 3;
offset_figlab = 1.5;
figspace_x = 5;
figspace_y = 2;


% Figure Eye data for revision of Seillier et al, 2017 =========================

% A: FF peak - blank
ax_new = axes(axPars, 'position', [xbegin+1 ybegin sq sq]);
% open figures
fig = openfig(['Z:\Corinna\SharedFigs\individuals\ff_peak-blank_4.fig'],'invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
xlim([-20 25])
ylim([-20 25])
set(gca,'XTick',[-20 0 20])
set(gca,'YTick', [-20 0 20])
message = sprintf('drug');
% xlabel(gca,'time from stimulus onset (ms)')
ylabel(gca, message)
message = sprintf(' Fano Factor \n (peak - blank)');
title(gca,message)
xlabel(gca,'baseline')
% offset_axis(0.05, axPars)
% set(gca,'XColor','w')

% axes(axPars,'position',[xbegin-offset_figlab ybegin+1.5*offset_figlab 1 1])
% title('a','fontsize',15)
% axis off



% legend
ax_new = axes(axPars, 'position', [xbegin-offset_figlab+0.5 ybegin-1.2*sq 2.5*sq 1*sq]);
fig = openfig(['Z:\Corinna\SharedFigs\individuals\legend_ff_peak-blank.fig'],'invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
axis off


 
savefig(['Z:\Corinna\SharedFigs\Figure_FFdiff.fig'])
print(h,'-dpdf',['Z:\Corinna\SharedFigs\Figure_FFdiff.pdf'],sprintf('-r%d',800))
% printeps(1, 'Z:\Corinna\SharedFigs\individuals\Figure_EyeData')

