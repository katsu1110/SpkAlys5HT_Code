function printfigureEye
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
ybegin = 16;
sq = 3;
offset_figlab = 1.5;
figspace_x = 5;
figspace_y = 2;


% Figure Eye data for revision of Seillier et al, 2017 =========================

% A: fixation precision
ax_new = axes(axPars, 'position', [xbegin ybegin-sq-figspace_y sq sq]);
% open figures
fig = openfig('Z:\Corinna\SharedFigs\fixationspan.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
xlim([0.2 1])
ylim([0.2 1])
set(gca,'XTick',[0.2 1])
set(gca,'YTick', [0.2 1])
message = sprintf('drug');
% xlabel(gca,'time from stimulus onset (ms)')
ylabel(gca, message)
message = sprintf(' fixation \n precision \n (dva^2)');
title(gca,message)
% offset_axis(0.05, axPars)
% set(gca,'XColor','w')

axes(axPars,'position',[xbegin-offset_figlab ybegin-figspace_y-1 1 1])
title('a','fontsize',15)
axis off

% B2: microsaccade amplitude
ax_new = axes(axPars, 'position', [xbegin+2*figspace_x ybegin-sq-figspace_y sq sq]);
% open figures
fig = openfig('Z:\Corinna\SharedFigs\microsac_amplitude.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
xlim( [0.03 0.12]);
ylim([0.03 0.12]);
set(gca,'XTick',[0.05 0.1])
set(gca,'YTick', [0.05 0.1])
message = sprintf('  microsaccade \n amplitude \n (dva)');
title(gca,message)
% message = sprintf('drug');
% xlabel(gca,'time from stimulus onset (ms)')
% ylabel(gca, message)
% ylabel(gca,'disparity ({^o})')
% set(gca, 'YTick', [-0.3 0 0.3])
% [a1,a2] = offset_axis(0.05, axPars);
% set(a1,'YTick',[])
% set(a1,'YColor','w')
% set(gca, 'XTick', [])
% set(gca,'XColor','w')


axes(axPars,'position',[xbegin+2*figspace_x-offset_figlab ybegin-figspace_y-1 1 1])
title('c','fontsize',15)
axis off

% B1: microsaccade rate
ax_new = axes(axPars, 'position', [xbegin+figspace_x ybegin-sq-figspace_y sq sq]);
% open figures
fig = openfig('Z:\Corinna\SharedFigs\microsac_counts.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
% caxis(crange)
% load('Z:\Katsuhisa\pupil_project\Figures\Figure10_kernelAvrew\individualFigs\avrew_rew1_high_tr.mat')
% text(0.95, 0.4, 'large', 'color', map(10,:),'fontsize',8, 'fontweight', 'bold')
% text(2, 0.4, strcat(' (n = ', num2str(tr2),')'), 'fontsize',8)
set(gca,'xscale','log')
set(gca,'yscale','log') 
xl = [0.5 30];
yl = [0.5 30];
xlim( xl );
ylim( yl );
set(gca,'XTick',xl, 'XTickLabel',{'0.3','30'})
set(gca,'YTick', yl, 'YTickLabel',{'0.3','30'})
xlabel(gca,'baseline')
message = sprintf('  microsaccade \n frequency');
title(gca,message)
% ylabel(gca,'disparity ({^o})')
% set(gca, 'XTick', [0.5 2.5 4.5],'XTickLabel',{'0','750','1500'})
% set(gca, 'YTick', [-0.3 0 0.3])
% 
% c = colorbar('eastoutside');
% cpos = c.Position;
% cpos(1) = 1.28*cpos(1);
% cpos(3) = 0.85*cpos(3);
% c.Position = cpos;
% c.AxisLocation = 'in';
% c.Box = 'off';
% c.FontSize = 8;
% c.TickDirection = 'out';
% c.Ticks = [-0.4 0 0.4];
% 
% offset_axis(0.05, axPars)
% set(gca,'XTick',[])
% set(gca,'XColor','w')

axes(axPars,'position',[xbegin+figspace_x-offset_figlab ybegin-figspace_y-1 1 1])
title('b','fontsize',15)
axis off

% legend
ax_new = axes(axPars, 'position', [xbegin+figspace_x-offset_figlab-0.5 ybegin-2.2*sq-figspace_y 2.5*sq 1*sq]);
fig = openfig('Z:\Corinna\SharedFigs\legend.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
axis off
% text(0.95, 0.4, 'large', 'color', map(10,:),'fontsize',8, 'fontweight', 'bold')
% text(2, 0.4, strcat(' (n = ', num2str(tr2),')'), 'fontsize',8)


 
savefig('Z:\Corinna\SharedFigs\Figure_EyeData.fig')
print(h,'-dpdf','Z:\Corinna\SharedFigs\Figure_EyeData.pdf',sprintf('-r%d',800))
% printeps(1, 'Z:\Corinna\SharedFigs\Figure_EyeData')

