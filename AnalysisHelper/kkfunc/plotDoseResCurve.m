function plotDoseResCurve

% open fig
close all
uiopen('Z:\Corinna\SerotoninPaper\Figures\Figure01_EffectMeanFiringRate\raw_figs\Recovery\DoseResponseFct.fig',1)

% get data from the fig
h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

% possible x vals
unix = [5 10 20 30 40];

% get y data in each x
for i = 1:5
    for k = 1:94
        if dataObjs(k).XData==unix(i)
            dose.x(i).y(k) = dataObjs(k).YData;
            if dataObjs(k).Marker == 'square'
                dose.x(i).animal(k) = 0;
            elseif dataObjs(k).Marker == 'o'
                dose.x(i).animal(k) = 1;
            end
        else
            dose.x(i).y(k) = nan;
            dose.x(i).animal(k) = nan;
        end
    end
    % remove nan
    dose.x(i).y = dose.x(i).y(~isnan(dose.x(i).y));
    dose.x(i).animal = dose.x(i).animal(~isnan(dose.x(i).animal));
    
    % recompute mean & SEM
    dose.x(i).mean = mean(exp(dose.x(i).y));
    dose.x(i).sem = std(exp(dose.x(i).y))/sqrt(length(dose.x(i).y));
end

figure;
for i = 1:5
    for k = 1:length(dose.x(i).y)
        if dose.x(i).animal(k)==0
            shape = 's';
        else
            shape = 'o';
        end
        scatter(unix(i), exp([dose.x(i).y(k)]), 20, 'filled', shape, ...
            'markerfacecolor', 'r', 'markeredgecolor', 'r', ...
            'markerfacealpha', 0.5, 'markeredgealpha', 0.8);
        hold on;
    end
end
hold on;
plot([unix(1) unix(5)], [1 1], '-k', 'linewidth', 0.5)
hold on;
errorbar(unix, [dose.x.mean], [dose.x.sem], 'color', [0 0.4 0], 'linewidth', 1)
hold on;
y_bottom = [dose.x.mean] - [dose.x.sem];
y_top = [dose.x.mean] + [dose.x.sem];
fill([unix fliplr(unix)],[y_bottom fliplr(y_top)], [0 0.4 0],'linewidth',0.5, ...
    'facecolor',[0 0.4 0], 'facealpha',0.1,...
    'edgecolor',[0 0 0], 'edgealpha',0.2);
 
xlim([5 40])
ylim([0.25 2])
set(gca,'YScale','log')
set(gca, 'YTick', [0.25 0.5 0.75 1 2])

