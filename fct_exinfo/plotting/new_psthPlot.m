function new_psthPlot( exinfo, exB, exD )
%plots a smoothed psth for every window
% NOTE: LATENCY IS OVERWRITTEN TO FULL POPULATION RESUTLS

h = figure('Name', exinfo.figname, 'UserData', exinfo, ...
    'Position', [526   227   953   754]);

col= getCol(exinfo);

[parB, latB, pvalB, sB] = plotHelper(exinfo, exB, 0, 0, 0.5, 'baseline');
[parD, latD, pvalD, sD] = plotHelper(exinfo, exD, 9, 6, 0.01, exinfo.drugname);

for i = 1:5
    sB(i).Position = ...
        [sB(i).Position(1)-0.05 ...
        sB(i).Position(2)+0.03 ...
        sB(i).Position(3:4)];
    sD(i).Position = [sD(i).Position(1)-0.05 ...
        sD(i).Position(2)-0.03 ...
        sD(i).Position(3:4)];
end

%%% All Data Plots
%stimuli depending latency
s(2) = subplot(3, 3, 6); % stimuli vs latency

plot(parB, latB(:, 3), 'o-', 'Color', col, ...
    'MarkerFaceColor', col); hold on;
plot(parD, latD(:, 3), 'o--', 'Color', col);

xlabel(exinfo.param1);  ylabel('latency');
s(2).XTick = parB;      box off;
legend('baseline', exinfo.drugname, 'Location','northoutside','Orientation','horizontal');
xlim([min([parB parD]), max([parB parD])]);

%latency regression
s(3) = subplot(3, 3, 9);  % latency baseline vs latency 5ht/nacl

scatter(latB(ismember(parB, parD), 3), ...
    latD(ismember(parD, parB), 3), col, 'filled');
eqax; regl; unity;

xlabel('baseline');     ylabel(exinfo.drugname);
title('all latency');   axis square

%TC
s(1) = subplot(3,3,3);

idx = exinfo.ratepar < 1000;
plot(exinfo.ratepar(idx), exinfo.ratemn(idx), 'o-', 'Color', col, 'MarkerFaceColor',  col);
hold on;
text(exinfo.ratepar(idx), exinfo.ratemn(idx) , ...
    ntrialwstar(exB.trials_n(idx), pvalB(idx, 3)), 'FontSize', 11);

idx = exinfo.ratepar_drug < 1000;
plot(exinfo.ratepar_drug(idx), exinfo.ratemn_drug(idx), 'o--', 'Color',  col);
text(exinfo.ratepar_drug(idx), exinfo.ratemn_drug(idx), ...
    ntrialwstar(exB.trials_n(idx), pvalD(idx, 3)), 'FontSize', 11);

xlim([min([parB parD]), max([parB parD])]);
xlabel(exinfo.param1); ylabel('spk rate');
s(1).XTick = exinfo.ratepar(idx);

legend('baseline', exinfo.drugname, 'Location','northoutside','Orientation','horizontal');
box off;

if strcmp(exinfo.param1, 'co') || strcmp(exinfo.param1, 'sf')
    set(s(1:2), 'XScale', 'log');
end


%% finals
savefig(h, exinfo.fig_psth);
close(h);
end


%%
function [parvls, lat, pval, s] = plotHelper(exinfo, ex, off, off2, off3, titletext)

s(1) = subplot(6, 3, 1+ off); % window 1
s(2) = subplot(6, 3, 4+ off); % window > 1
s(3) = subplot(6, 3, 7+ off); % all


s(4) = subplot(4, 3, 2+ off2); % param1 vs latency
s(5) = subplot(4, 3, 5+ off2); % window 1 vs  window > 1

allTrials = addWindow(exinfo, ex.Trials);

% stimuli parameter
parvls = unique([allTrials.(exinfo.param1)]);
parvls = parvls(parvls < 1000);

% predefine variables
col = lines(length(parvls));
lat = nan(length(parvls), 3);
n = nan(length(parvls), 2);
pval = lat;

% loop through all stimuli parameter and get psth and latency
for pari = 1:length(parvls)
    
    ind1 = [allTrials.(exinfo.param1)] == parvls(pari);
    ind2 = [allTrials.window] == 1;
    
    % window == 1
    [lat(pari, 1), pval(pari, 1), psth1, n(pari, 1)] = getLatencyDG(exinfo, allTrials(ind1&ind2), true);
    
    % window > 1
    [lat(pari, 2), pval(pari, 2), psth2, n(pari, 2)] = getLatencyDG(exinfo, allTrials(ind1&~ind2), true);
    
    % all
    [lat(pari, 3), pval(pari, 3), psth_all, ~] = getLatencyDG(exinfo, allTrials(ind1), true);
    
    % plotting
    axes(s(1))
    plot(psth1, 'Color', col(pari, :));  hold on;
    axes(s(2))
    plot(psth2, 'Color', col(pari, :)); hold on
    axes(s(3))
    plot(psth_all, 'Color', col(pari, :), 'LineWidth', 1.5); hold on
    
end

% plot latency across stimuli parameter
axes(s(4))
plot(parvls, lat(:, 3), 'x-', 'Color', [0 0 0 0.5]); hold on;
plot(parvls, lat(:, 1), 'bo-', 'MarkerFaceColor', 'b'); hold on;
plot(parvls, lat(:, 2), 'bo--'); hold on;

legend('all', 'seq# == 1', 'seq# > 1', 'Location','northoutside','Orientation','horizontal');
xlabel(exinfo.param1); ylabel('latency');
s(4).XTick = parvls; s(4).XLim(2) = max(parvls);

for pari = 1:length(parvls)
    if pval(pari, 1) < 0.05 && pval(pari, 2) < 0.05
        text(parvls(pari),s(4).YLim(2), sprintf('%1i */%1i *', n(pari, :)), ...
            'FontSize', 8, 'FontWeight', 'bold');
    elseif pval(pari, 1) < 0.05
        text(parvls(pari),s(4).YLim(2), sprintf('%1i */%1i', n(pari, :)), ...
            'FontSize', 8, 'FontWeight', 'bold');
    elseif pval(pari, 1) < 0.05
        text(parvls(pari),s(4).YLim(2), sprintf('%1i/%1i *', n(pari, :)), ...
            'FontSize', 8, 'FontWeight', 'bold');
    else
        text(parvls(pari),s(4).YLim(2), sprintf('%1i/%1i', n(pari, :)), ...
            'FontSize', 8, 'FontWeight', 'bold');
    end
end
box off;


if strcmp(exinfo.param1, 'co') || strcmp(exinfo.param1, 'sf')
    set(s(4), 'XScale', 'log');
end

%plot latency wind=1 vs latency wind>1
axes(s(5))
scatter(lat(:, 1), lat(:, 2), 'o', 'filled');
eqax; regl; unity;
axis square
xlabel('seq# == 1');
ylabel('seq# > 1');
s(5).Position = [s(5).Position(1:2) s(5).Position(3:4)-0.02];

%%%
axes(s(1))
rightAdjTitle('seq# == 1');
s(1).Position(4) = s(1).Position(4) - 0.02;


axes(s(2))
rightAdjTitle('seq# > 1');
ylabel(titletext, 'FontSize', 12, 'Fontweight', 'bold');
s(2).Position(4) = s(2).Position(4) - 0.02;


axes(s(3))
rightAdjTitle('all');
xlabel('time');
s(3).Position(4) = s(3).Position(4) - 0.02;



%%%
set(s(1:3), 'YLim', ...
    [ min( min(vertcat(s(1:3).YLim)) ), max( max(vertcat(s(1:3).YLim)) )], ...
    'XLim', [0 450], 'Box', 'off');
% set(s(1:2), 'XTick', []);
set(s(1:5), 'TickDir', 'out');


% parameters
axes('Position', [0.02 off3 0.5 0.01]);
l = length(parvls);
xlim([0, l+1]);

for pos = 1:l
    plot([pos-0.5, pos+0.5], [0 0], 'Color', [col(pos, :) 0.5], 'LineWidth', 3);
    hold on;
    text(pos-0.35, 0, sprintf('%1.2f ', parvls(pos)), ...
        'FontSize', 6);
end
axis off
end



function rightAdjTitle(titletext)

t = title(titletext);

set(t, 'horizontalAlignment', 'right');

set(t, 'units', 'normalized');

h1 = get(t, 'position');

set(t, 'position', [1 h1(2) h1(3)]);

end



function c = ntrialwstar(n , p)

c = cell(length(n), 1);

for i  =1:length(n)
    
    c{i} = num2str(n(i));
    if p(i) <0.5
        c{i} = [c{i} '*'];
    end
end

end

