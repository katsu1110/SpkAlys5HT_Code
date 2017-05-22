function h = znormplot( ex_base, ex_drug, exinfo )
% znormplot
% plots the z-scored data in the drug and baseline condition 


h = figure('Name', exinfo.figname);
drugcolor = getCol(exinfo);

% 1. baseline condition
ax(1) = subplot(9,1,[3 4 5]);
s1 = plotZScoreHelper(ex_base, true);% baseline
s1.Color = drugcolor; crossl;


title(sprintf(['Fano Factor mean all and for 2nd half \n'...
    ' Baseline: %1.2f, %1.2f \n '],...
    nanmean(exinfo.ff.classic.ff),  nanmean(exinfo.ff.classic_2ndhalf.ff)), ...
    'FontSize', 9); 

% 2. drug condition
ax(2) = subplot(9,1,[7 8 9]);
s2 = plotZScoreHelper(ex_drug, false);% drug
s2.Color = drugcolor; s2.LineStyle = ':';crossl;

title(sprintf(['Fano Factor mean all and for 2nd half \n'...
    exinfo.drugname ': %1.2f, %1.2f '],...
    nanmean(exinfo.ff_drug.classic.ff),  nanmean(exinfo.ff_drug.classic_2ndhalf.ff)),...
    'FontSize', 9); 


% 3. set labels and limits of y-axes equal
xlabel(ax(2), '# stimulus presentation');
ylabel(ax(1), 'z-normed spike rate');


allYLim = get(ax, {'YLim'});
allYLim = cat(2, allYLim{:});
set(ax, 'YLim', [min(allYLim), max(allYLim)]);



% 4. plot stimulus color coding
axes('position', [0.1 0.9 0.9 0.1]);
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




% 5. save figure
savefig(h, exinfo.fig_varxtime);
delete(h);

end



%%
function s = plotZScoreHelper(ex, fill_flag)
% plot z-normed spike rate over time

% plot the trend line
x = 1:length([ex.Trials]);
y = [ex.Trials.zspkcount];
s = plot(x, y, 'k-', 'LineWidth', 0.5); hold on;

% plot trials with equal stimulus in the same color
[ stimparam, vals] = getStimParam( ex );
nvals = length(vals);
col = lines(nvals);

for i = 1:nvals
    ind = [ex.Trials.(stimparam)] == vals(i);
    col_i = col(i,:);
    if fill_flag
        scatter(x(ind), y(ind), 20, col_i, 'filled'); hold on;
    else
        scatter(x(ind), y(ind), 20, col_i); hold on;
    end
end
end