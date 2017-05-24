function exinfo = new_psthPlot_red( exinfo, exB, exD )
%plots a smoothed psth for every window
% NOTE: LATENCY IS OVERWRITTEN TO FULL POPULATION RESUTLS

h = figure('Name', exinfo.figname, 'UserData', exinfo, ...
    'Position', [897   588   682   396]);

col= getCol(exinfo);

[parB, latB, pvalB, sB, ntrialsB] = plotHelper(exinfo, exB, 0, 'baseline');
[parD, latD, pvalD, sD, ntrialsD] = plotHelper(exinfo, exD, 3, exinfo.drugname);

sB(1).Position(1) = 0.08;
sB(1).Position(2) = sB(1).Position(2)+0.05;
sB(1).Position(4) = sB(1).Position(4)-0.1;

sD(1).Position(1) = 0.08;
sD(1).Position(2) = sD(1).Position(2)+0.05;
sD(1).Position(4) = sD(1).Position(4)-0.1;

%%% All Data Plots
%stimuli depending latency
s = subplot(1, 3, [2 3]); % stimuli vs latency

plot(parB, latB, 'o-', 'Color', col, 'MarkerFaceColor', col); hold on;
plot(parD, latD, 'o--', 'Color', col);

xlabel(exinfo.param1);  ylabel('latency');
s.XTick = parB;      box off;
legend('baseline', exinfo.drugname, 'Location','northoutside','Orientation','horizontal');
xlim([min([parB parD]), max([parB parD])]);


text(parB, latB, ntrialwstar(ntrialsB, pvalB), 'FontSize', 11);
text(parD, latD, ntrialwstar(ntrialsD, pvalD), 'FontSize', 11);

if strcmp(exinfo.param1, 'co') || strcmp(exinfo.param1, 'sf')
    set(s, 'XScale', 'log');
end


%% finals
savefig(h, exinfo.fig_psth);
close(h);


%% add latency to exinfo

exinfo.lat = [parB; latB'; pvalB'];
exinfo.lat_drug = [parD; latD'; pvalD'];

end



%%
function [parvls, lat, pval, s, ntrials] = plotHelper(exinfo, ex, off, titletxt)

s(1) = subplot(2, 3, 1+off); % all

allTrials = addWindow(exinfo, ex.Trials);

% stimuli parameter
parvls = unique([allTrials.(exinfo.param1)]);
parvls = parvls(parvls < 1000);

% predefine variables
col = lines(length(parvls));
lat = nan(length(parvls), 1);
pval = lat;

% loop through all stimuli parameter and get psth and latency
for pari = 1:length(parvls)
    
    ind = [allTrials.(exinfo.param1)] == parvls(pari);
    
    % all
    [lat(pari), pval(pari), psth_all, ~] = getLatencyDG(exinfo, allTrials(ind), true);
    
    % plotting
    plot(psth_all, 'Color', col(pari, :), 'LineWidth', 1.5); hold on
    
    ntrials(pari) = sum(ind);
end

box off;

axes(s(1))
rightAdjTitle(['all ' titletxt]); 
xlabel('time');
ylabel('spk/s');


%%%
set(s, 'XLim', [0 450], 'Box', 'off', 'TickDir', 'out');

% parameters
s(2) = axes('Position', [0.01 s(1).Position(2)-0.1 0.3 0.03]);
l = length(parvls);
xlim([0, l+1]);

for pos = 1:l
    plot([pos-0.5, pos+0.5], [0 0], 'Color', [col(pos, :) 0.5], 'LineWidth', 3);
    hold on;
    text(pos-0.35, 0, sprintf('%1.2f ', parvls(pos)), ...
        'FontSize', 8);
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
    if p(i) <0.001
        c{i} = [c{i} '***'];
    elseif p(i) <0.01
        c{i} = [c{i} '**'];
    elseif p(i) <0.5
        c{i} = [c{i} '*'];
    end
end

end

