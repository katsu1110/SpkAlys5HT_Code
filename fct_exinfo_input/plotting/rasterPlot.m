function rasterPlot( exinfo, ex, ex_drug )
% classical raster plot for both conditions with stimulus indicating color 
% and lines indicating trial that were first in the fixation period
% 
% saves the figure in exinfo.fig_raster
% 
% @CL 16.11.2016
% added the development of spike rate over time as additional axes
% @CL 06.02.2017


% --------------------------------------- plot
h = figure('Name', exinfo.figname, 'UserData', exinfo, ...
    'Position', [680   156   560   822]);
rasterPlotHelper(exinfo, ex, ex_drug)

% --------------------------------------- save
savefig(h, exinfo.fig_raster);
close(h);
end



%%
function rasterPlotHelper(exinfo, ex, ex_drug)
% this is the meta function for the true plotting
% it calls to plot the raster for each condition with plotCondition(..) 
% and the stimulus annotation in the lower part of the figure

if exinfo.isadapt
    time = 0:0.001:5;
else
    time = -0.05:0.001:0.45;
end

% 1. plot baseline condition
ntrials = length([ex.Trials])+1;

% baseline raster
ax_rast = subplot(9,7,[(8:7:63), (9:7:63)]);
plotCondition(addWindow(exinfo, ex.Trials), time, exinfo.param1, exinfo.ratepar)
title('baseline'); ylim([1 ntrials]); ax_rast.View =[0 -90];

% baseline rate development
ax = subplot(9,7,(10:7:63));
plotRateHelper(ex, exinfo.ff, exinfo.ratepar);
ax.View = [90 90];  ax.XTick = []; 
ylabel('spk/s');    xlim([1 ntrials]);
box off

% 2. plot drug condition
ntrials = length([ex_drug.Trials])+1;
% drug raster
ax_rast = subplot(9,7,[(12:7:63), (13:7:63)]);
plotCondition(addWindow(exinfo, ex_drug.Trials), time, exinfo.param1, exinfo.ratepar)
title(exinfo.drugname); ylim([1 ntrials]);ax_rast.View =[0 -90];

% drug  rate development
ax = subplot(9,7,(14:7:63));
plotRateHelper(ex_drug, exinfo.ff_drug, exinfo.ratepar_drug);
ax.View = [90 90];  ax.XTick = [];      
ylabel('spk/s');    xlim([1 ntrials]);
box off


% 3. plot stimulus parameter and corresponding color
subplot(9,7,1:7);
l = length(exinfo.ratepar); 
xlim([0, l+1]); 
col = lines(l);
for pos = 1:l
    plot([pos-0.5, pos+0.5], [0 0], 'Color', [col(pos, :) 0.5], 'LineWidth', 3);
    hold on;
    text(pos-0.35, 0, sprintf('%1.2f ', exinfo.ratepar(pos)), 'FontSize', 9);
end
text(0.5, 0.5, [exinfo.param1 ', grey=1st trial'], 'FontSize', 9);
axis off;
end



function plotCondition(Trials, time, param, parvls)
%%% plot the raster matrix according to different stimuli
% also indicate phase and window 

row_i = 1;
col = lines(length(parvls));

% all the trials 
for par_i = 1:length(parvls)
    trials = Trials( [Trials.(param)] == parvls(par_i));
    
    row_i = getRasterPlotMatrix(trials, time, row_i, col(par_i, :));
    plot(time([1, end]), [row_i, row_i], 'Color', [.5 .5 .5]);
    
end

xlim([-0.05 0.45]);
plot([0 0], [1, row_i], 'Color', [0.5 0.5 0.5]);

set(gca, 'Clipping', 'off', 'TickDir', 'out'); box off;
end


function [y_idx, raster] = getRasterPlotMatrix(Trials, time, y_idx, col)
% converts the time stemps of spikes into a sequence of 0s and 1s, with 1s
% indicating a spike event. Each row is one trial.

% initiate raster
raster = nan(length(time), length(Trials));

% each row is a trial, each column a time bin
for t = 1:length(Trials)
    
    t_strt = Trials(t).Start - Trials(t).TrialStart;
    
    if time(end)>0.5
        t_strt = t_strt(t_strt<=5);
    end
    
    % spikes in the time range
    spk =  Trials(t).Spikes(...
        Trials(t).Spikes>=t_strt(1)+min(time) &  Trials(t).Spikes<=t_strt(end));
    spk = spk-t_strt(1); % align spikes to stimulus onset
    
    % adding each spike event as '1' to the raster matrix
    idx = round(( spk + abs(time(1)) ) *1000);
    idx(idx==0) = 1; % avoid an index of 0
    raster(idx, t) = 1;
    
    % plotting
    if Trials(t).window == 1
        fill([time(1), time(end), time(end), time(1)],...
            [y_idx, y_idx, y_idx+0.95, y_idx+0.95], ...
            [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeAlpha', 0); hold on;
    end
    plot([spk spk]', repmat([y_idx; y_idx+1], 1, length(idx)), ...
        'Color', col);
    plot([time(end); time(end)], [y_idx; y_idx+0.95], 'Color', col);
    
    y_idx = y_idx+1;
    hold on;

end

end



function plotRateHelper(ex, ff, ratepar)

x = 1.5;
y = [ex.Trials.spkRate];

% plot each stimulus response individually with stimulus specific color
% plot trials with equal stimulus in the same color
[ stimparam, vals] = getStimParam( ex );
nvals = length(vals);
col = lines(nvals);

for i = 1:nvals
    ind = [ex.Trials.(stimparam)] == vals(i);
    col_i = col(i,:);
    plot(x:x+sum(ind)-1, y(ind), '-', 'Color', col_i); hold on;
    scatter(x:x+sum(ind)-1, y(ind), 15, col_i, 'filled'); hold on;
    idx = find(ratepar == vals(i));
    
    try
        text(x+sum(ind)/2, max(y), ...
            sprintf('ff(%1.0f)= %1.1f / %1.1f\n=%1.2f \n ff2(%1.0f)= %1.1f / %1.1f\n=%1.2f',...
            ff.classic.stimrep(idx), ....
            ff.classic.spkcnt_var(idx), ff.classic.spkcnt_mn(idx), ff.classic.ff(idx),...
            ff.classic_2ndhalf.stimrep(idx), ....
            ff.classic_2ndhalf.spkcnt_var(idx), ff.classic_2ndhalf.spkcnt_mn(idx), ...
            ff.classic_2ndhalf.ff(idx)), 'FontSize', 7, 'Color', col_i);
        
    catch
        text(x+sum(ind)/2, max(y), ...
            sprintf('ff= %1.1f / %1.1f\n=%1.2f',...
            ff.classic.spkcnt_var(idx), ff.classic.spkcnt_mn(idx), ff.classic.ff(idx)),...
            'FontSize', 7, 'Color', col_i);
    end
    x=x+sum(ind);
    plot([x x]-0.5, [0, 150], 'Color', [.5 .5 .5]);
end

ylim([0 max(y)]);

end