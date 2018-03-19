function h = plotIntr(exinfo)
% h is the figure handle of the figure containing the superimposed raw and
% fitted tuning curves for both experiments +/- sme using the information
% from exinfo.fitparam and exinfo.fitparam_drug. the data is provided by
% evalSingleDG where, depending on the stimulus feature, descriptive
% functions were fitted to the raw tuning curves.
% 
% the figure is saved as exinfo.fig_tc
% 
% written by Katsuhisa (07.06.17)
%++++++++++++++++++++++++++++++++++++++++++++

h = figure('Name', exinfo.figname);

% tuning curve: no drug
subplot(2,6,1)
xv1 = exinfo.ratepar_1(1:end-1);
yv1 = exinfo.ratemn_1(1:end-1);
ye1 = exinfo.ratesme_1(1:end-1);
xv2 = exinfo.ratepar_2(1:end-1);
yv2 = exinfo.ratemn_2(1:end-1);
ye2 = exinfo.ratesme_2(1:end-1);
niceerrorbar(xv1, yv1, ye1, {'-b'}); hold on;
% plot(xfit1, yfit1, '-b'); hold on;
niceerrorbar(xv2, yv2, ye2, {'-r'}); hold on;
title('baseline')
hold on;
p1 = plot(xv1, yv1, '-b');
hold on;
p2 = plot(xv2, yv2, '-r');
legend([p1,p2],'small ps','large ps','location','best')
legend('boxoff')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
yy0 = get(gca, 'YLim');

    % plot(xfit2, yfit2, '-r'); hold on;
% axis square
    
% tuning curve: drug
subplot(2,6,2)
xv1 = exinfo.ratepar_1_drug(1:end-1);
yv1 = exinfo.ratemn_1_drug(1:end-1);
ye1 = exinfo.ratesme_1_drug(1:end-1);
xv2 = exinfo.ratepar_2_drug(1:end-1);
yv2 = exinfo.ratemn_2_drug(1:end-1);
ye2 = exinfo.ratesme_2_drug(1:end-1);
niceerrorbar(xv1, yv1, ye1, {'-','color','b'}); hold on;
% plot(xfit1, yfit1, '-b'); hold on;
niceerrorbar(xv2, yv2, ye2, {'-','color','r'}); hold on;
switch exinfo.is5HT
    case 1
        drugname = '5HT';
        col = [1 0 0];
    otherwise
        drugname = 'NaCl';
        col = [0 0 0];
end
title(drugname)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
yy2 = get(gca, 'YLim');
% plot(xfit2, yfit2, '-r'); hold on;

subplot(2,6,1)
ylim([min([yy0 yy2]), max([yy0, yy2])])
subplot(2,6,2)
ylim([min([yy0 yy2]), max([yy0, yy2])])

% glm weights
subplot(2,6,3)
% bar(1:3, [exinfo.b_pf_drug exinfo.b_upf_drug;...
%         exinfo.b_pf_pupil exinfo.b_upf_pupil;...
%         exinfo.b_pf_intr exinfo.b_upf_intr]);
bar(1:3, exinfo.psglm.b(3, 1:3),...
    'FaceColor',col,'EdgeColor','w')
set(gca, 'XTick', 1:3, 'XTickLabel', {'drug','pupil','interaction'})
xtickangle(45)
xlim([0.5 3.5])
ylabel('beta weight (a.u.)')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% legend('prefered','unprefered','location','eastoutside')
% legend('boxoff')

% glm model performance
subplot(2,6,4)
bar(1:4, exinfo.psglm.perf(1:4),...
    'FaceColor',col,'EdgeColor','w')
set(gca, 'XTick', 1:4, 'XTickLabel', {'stm','+drug','+pupil','+interaction'})
xtickangle(45)
xlim([0.5 4.5])
ylabel('model performance (a.u.)')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% pupil size time-course
subplot(2,6,5)
nstm = unique(exinfo.psglm.tc(:,2));
nlen = length(nstm);
if nlen==1
    offset = 0;
else
    offset = 0.1;
end
ncol = size(exinfo.psglm.tc,2)-2;
for d = 1:2
    switch d
        case 1
            colt = 'b';
        case 2
            colt = col;
    end
    for n = 1:length(nstm)
        v = exinfo.psglm.tc(exinfo.psglm.tc(:,2)==nstm(n)...
            & exinfo.psglm.tc(:,1)==d-1, 3:end);

        plot(([[1:ncol] + ncol*(n-1)]/500) + offset*(n-1), mean(v,1), '-', 'color', colt) 
        hold on;
    end
end
xlabel('time (s)')
ylabel('pupil size (a.u.)')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% interaction table
subplot(2,6,6)
imagesc(exinfo.psglm.inter_table)
colormap(jet);
set(gca, 'XTick', 1:2, 'XTickLabel', {'5HT','baseline'})
set(gca, 'YTick', 1:2, 'YTickLabel', {'S-ps','L-ps'})
text(0.7,1,num2str(exinfo.psglm.inter_table(1,1)),'color','w')
text(0.7,2,num2str(exinfo.psglm.inter_table(2,1)),'color','w')
text(1.7,2,num2str(exinfo.psglm.inter_table(2,2)),'color','w')
text(1.7,1,num2str(exinfo.psglm.inter_table(1,2)),'color','w')

if ~exinfo.isRC
    % p-anova
    subplot(2,6,7)
    bar([exinfo.p_anova_1 exinfo.p_anova_2 exinfo.p_anova_1_drug exinfo.p_anova_2_drug],...
        'FaceColor',col,'EdgeColor','w')
    set(gca, 'XTick', 1:4, 'XTickLabel', {'smallPSxBase','largePSxBase','smallPSxDrug','largePSxDrug'})
    xtickangle(45);
    title('p-anova')
    xlim([0.5 4.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    % ff
    subplot(2,6,8)
    bar([nanmedian(exinfo.ff_1.classic.ff) nanmedian(exinfo.ff_2.classic.ff) ...
           nanmedian(exinfo.ff_1_drug.classic.ff) nanmedian(exinfo.ff_2_drug.classic.ff)],...
        'FaceColor',col,'EdgeColor','w')
    title('ff')
    set(gca, 'XTick', 1:4, 'XTickLabel', {})
    xlim([0.5 4.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    % legend('ps-small, no drug','ps-large, no drug','ps-small, drug','ps-large, drug', 'location','best')
    % legend('boxoff')

    % nc
    subplot(2,6,9)
    bar([exinfo.rsc_1 exinfo.rsc_2 exinfo.rsc_1_drug exinfo.rsc_2_drug],...
        'FaceColor',col,'EdgeColor','w')
    title('nc')
    set(gca, 'XTick', 1:4, 'XTickLabel', {})
    xlim([0.5 4.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    % legend('ps-small, no drug','ps-large, no drug','ps-small, drug','ps-large, drug', 'location','best')
    % legend('boxoff')

    % sc
    subplot(2,6,10)
    bar([exinfo.rsig_1 exinfo.rsig_2 exinfo.rsig_1_drug exinfo.rsig_2_drug],...
        'FaceColor',col,'EdgeColor','w')
    title('sc')
    set(gca, 'XTick', 1:4, 'XTickLabel', {})
    xlim([0.5 4.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    % legend('ps-small, no drug','ps-large, no drug','ps-small, drug','ps-large, drug', 'location','best')
    % legend('boxoff')

    % tuning curve height
    subplot(2,6,11)
    bar([exinfo.tcdiff_1 exinfo.tcdiff_2 exinfo.tcdiff_1_drug exinfo.tcdiff_2_drug],...
        'FaceColor',col,'EdgeColor','w')
    title('tc_height')
    set(gca, 'XTick', 1:4, 'XTickLabel', {})
    xlim([0.5 4.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    % legend('ps-small, no drug','ps-large, no drug','ps-small, drug','ps-large, drug', 'location','best')
    % legend('boxoff')

    % phase selectivity
    subplot(2,6,12)
    bar([exinfo.phasesel_1 exinfo.phasesel_2 exinfo.phasesel_1_drug exinfo.phasesel_2_drug],...
        'FaceColor',col,'EdgeColor','w')
    title('phase_sel')
    set(gca, 'XTick', 1:4, 'XTickLabel', {})
    xlim([0.5 4.5])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    % legend('ps-small, no drug','ps-large, no drug','ps-small, drug','ps-large, drug', 'location','best')
    % legend('boxoff')
end
       
set(gcf, 'position',[344         480        1451         498])
set(h, 'UserData', exinfo);
savefig(h, exinfo.fig_intr);
close(h);

end

