function  psthPlot( exinfo, exB, exD )
%plots a smoothed psth for every window
% NOTE: LATENCY IS OVERWRITTEN TO FULL POPULATION RESUTLS


h = figure('Name', exinfo.figname, 'UserData', exinfo, ...
    'Position', [526   227   953   754]);

% plots per window
s(1:5) = plotHelper(exinfo.psth, 'base', exinfo.ratepar(exinfo.ratepar<1000), ...
    0, exinfo.lat);
s(6:10) = plotHelper(exinfo.psth_drug, exinfo.drugname, ...
    exinfo.ratepar_drug(exinfo.ratepar_drug<1000), ...
    6, exinfo.lat_drug);

% averaged across all windows
[s(11:13) , latB, latD] = plotHelperAll(exinfo, exB, exD );

% assignment of averaged latencies 
exinfo.lat = latB;
exinfo.lat_drug = latD;

% align psth plots
ymax = max([s([1:4, 6:9 11:12]).YLim]);
set(s([1:4, 6:9 11:12]), 'ylim', [0 ymax]);

% align latency plots
ymax = max([s([5 10 13]).YLim]);
set(s([5 10 13]), 'ylim', [0 ymax]);


s(10).YLabel.String = '';
s(1).YLabel.String = '10ms';
s(7).YLabel.String = 'smoothed';

%use the correct parmaeter scale
if strcmp(exinfo.param1, 'co') || strcmp(exinfo.param1, 'sf')
    set(s([5 10 13]), 'XScale', 'log');
end

savefig(h, exinfo.fig_psth);
close(h);
end


%%
function s = plotHelper(psth, name, par, row_i, lat)

col = lines(length(par));

if row_i == 0
    s(5) = subplot(4, 6, [13 14 19 20] );
else
    s(5) = subplot(4, 6, [15 16 21 22] );
end
    

% smooth the psth with a gaussian kernel
kernel = gaussian(0, 20, 1, 0, -70:70); kernel = kernel/sum(kernel);

for wind_i = 1:4

    % smoothed psth
    s(wind_i) = subplot(4, 6, wind_i+row_i);
    
    for par_i = 1:length(par)
        plot(conv(psth{par_i, wind_i, 1}, kernel, 'same'), 'Color', col(par_i, :));
        hold on;
    end
    hold off
    xlim([0 450]);
    title(sprintf(['w %1d ' name ' psth'], wind_i));
    
    % latency and number of rep
    axes(s(5))
    plot(par, lat(:, wind_i), 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', wind_i);hold on;
    for par_i = 1:length(par)
         plot(par(par_i), lat(par_i, wind_i), 'o', ...
             'MarkerFaceColor', col(par_i, :), 'MarkerEdgeColor', col(par_i, :));
         hold on;
         text(par(par_i), lat(par_i, wind_i)+5, num2str(psth{par_i, wind_i, 2}), ...
         'FontSize', 10+wind_i/2, 'FontWeight', 'bold');
    end
    title([name ' thin-w1, thick-w4']);
    ylabel('latency');
    xlabel('stimulus parameter');
    set(gca, 'XTick', par);
    xlim([min(par) max(par)]);
end


end




function [s, latB, latD] = plotHelperAll(exinfo, exB, exD)


% smooth the psth with a gaussian kernel
kernel = gaussian(0, 20, 1, 0, -70:70); kernel = kernel/sum(kernel);
col = lines(max( [length(exinfo.ratepar) length(exinfo.ratepar_drug) ] ));


%%% smoothed psth for BASELINE
s(1) = subplot(4, 6, [5 6]);
paramB = exinfo.ratepar(exinfo.ratepar<1000);

for idx = 1:length(paramB)
    ctrials = exB.Trials ( [exB.Trials.(exinfo.param1)] == paramB(idx) );
    [latB(idx) psthB nB(idx)] = getLatencyDG( exinfo, ctrials, true );
    
    plot(psthB,  'Color', col(idx, :));
    hold on;
end
title('baseline');  xlim([0 450]);    
legend(num2str(paramB'), 'Location', 'EastOutside');



%%% smoothed psth for DRUG
s(2) = subplot(4, 6, [11 12]);
paramD = exinfo.ratepar_drug(exinfo.ratepar_drug<1000);

for idx = 1:length(paramD)

    ctrials = exD.Trials ( [exD.Trials.(exinfo.param1)] == paramB(idx) );
    [latD(idx) psthD nD(idx)] = getLatencyDG( exinfo, ctrials, true );
    
    plot(psthD, 'Color', col(idx, :));
    hold on;
end
title(exinfo.drugname); xlim([0 450]);
legend(num2str(paramD'), 'Location', 'EastOutside');



%%% latency and number of rep
s(3) = subplot(4, 6, [17 18 23 24]);
c = getCol(exinfo);

% Baseline
plot(paramB, latB, 'o-', 'Color', c, 'MarkerFaceColor', c);
for idx = 1:length(paramB)
    text(paramB(idx), latB(idx), num2str(nB(idx)), 'FontSize', 12);
    hold on;
end

% Drug
plot(paramD, latD, 'o--', 'Color', c);
for idx = 1:length(paramD)
    text(paramD(idx), latD(idx), num2str(nD(idx)), 'FontSize', 12);
    hold on;
end

title('baseline-, drug--');
xlabel(exinfo.param1);
set(gca, 'XTick', exinfo.ratepar, 'XLim', [min([paramB, paramD]) max([paramB, paramD])]); 

end



