function printfigureNCFF(minss)
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


% Figure 9 for revision of Seillier et al, 2017 =========================

% A: noise correlation
ax_new = axes(axPars, 'position', [xbegin+1 ybegin sq sq]);
% open figures
fig = openfig(['Z:\Corinna\SharedFigs\individuals\nc_' num2str(minss) '.fig'],'invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
xlim([-1 1])
ylim([-1 1])
set(gca,'XTick',[-1 0 1])
set(gca,'YTick', [-1 0 1])
message = sprintf('drug');
% xlabel(gca,'time from stimulus onset (ms)')
ylabel(gca, message)
message = sprintf(' Noise \n Correlation');
title(gca,message)
xlabel(gca,'baseline')
% offset_axis(0.05, axPars)
% set(gca,'XColor','w')

axes(axPars,'position',[xbegin-offset_figlab ybegin+1.5*offset_figlab 1 1])
title('a','fontsize',15)
axis off

% % B2: microsaccade amplitude
% ax_new = axes(axPars, 'position', [xbegin+2*figspace_x ybegin-sq-figspace_y sq sq]);
% % open figures
% fig = openfig('Z:\Corinna\SharedFigs\individuals\microsac_amplitude.fig','invisible');
% ax_old = findobj(fig, 'type', 'axes');
% copyobj(ax_old.Children, ax_new); delete(fig);
% xlim( [0.03 0.12]);
% ylim([0.03 0.12]);
% set(gca,'XTick',[0.05 0.1])
% set(gca,'YTick', [0.05 0.1])
% message = sprintf('  microsaccade \n amplitude \n (dva)');
% title(gca,message)
% % message = sprintf('drug');
% % xlabel(gca,'time from stimulus onset (ms)')
% % ylabel(gca, message)
% % ylabel(gca,'disparity ({^o})')
% % set(gca, 'YTick', [-0.3 0 0.3])
% % [a1,a2] = offset_axis(0.05, axPars);
% % set(a1,'YTick',[])
% % set(a1,'YColor','w')
% % set(gca, 'XTick', [])
% % set(gca,'XColor','w')


% axes(axPars,'position',[xbegin+2*figspace_x-offset_figlab ybegin-figspace_y-1 1 1])
% title('c','fontsize',15)
% axis off

% B: fano factor
ax_new = axes(axPars, 'position', [xbegin+1 ybegin-2*sq-figspace_y sq sq]);
% open figures
fig = openfig(['Z:\Corinna\SharedFigs\individuals\ff_' num2str(minss) '.fig'],'invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
xl = [0 15];
yl = [0 15];
xlim( xl );
ylim( yl );
set(gca,'XTick',xl, 'XTickLabel',{'0','15'})
set(gca,'YTick', yl, 'YTickLabel',{'0','15'})
message = sprintf(' Fano \n Factor');
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

axes(axPars,'position',[xbegin-offset_figlab ybegin-2*sq-figspace_y+1.5*offset_figlab 1 1])
title('b','fontsize',15)
axis off

% legend 1
ax_new = axes(axPars, 'position', [xbegin-offset_figlab ybegin-1.2*sq 2.5*sq 1*sq]);
fig = openfig(['Z:\Corinna\SharedFigs\individuals\legend_nc' num2str(minss) '.fig'],'invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
axis off

% legend 2
ax_new = axes(axPars, 'position', [xbegin-offset_figlab ybegin-2.5*sq-2*figspace_y 2.5*sq 1*sq]);
fig = openfig(['Z:\Corinna\SharedFigs\individuals\legend_ff' num2str(minss) '.fig'],'invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
axis off


 
savefig(['Z:\Corinna\SharedFigs\Figure_NCFF_' num2str(minss) '.fig'])
print(h,'-dpdf',['Z:\Corinna\SharedFigs\Figure_NCFF_' num2str(minss) '.pdf'],sprintf('-r%d',800))
% printeps(1, 'Z:\Corinna\SharedFigs\individuals\Figure_EyeData')

