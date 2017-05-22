function h = rcPlot( exinfo )
%specified plotting of sdfs when generating exinfo
%
%make two subplots for each baseline and drug condition
%
% @CL 22.01.2016


h = figure('Name', exinfo.figname, 'UserData', exinfo, 'Position', [680 274 560 704]);
rcplot_helper(exinfo)
set(findobj('type', 'axes'), 'fontsize', 8)

%--------------------------------- save
savefig(h, exinfo.fig_sdfs);
close(h);
end

%%
function rcplot_helper(exinfo)

c = hsv(length(exinfo.sdfs.s));
g = [0.9 0.9 0.9];
%------------------------------------------ baseline
baseplot = subplot(3,1,1);

% stimulus triggered SDFs 
for i = 1:length(exinfo.sdfs.s) % orientation
    plot(exinfo.times/10, exinfo.sdfs.s{i}, 'Color', c(i,:)); ho
    leg{i} = sprintf('%1.0f \t n=%1.0f', exinfo.sdfs.x(i), exinfo.sdfs.n(i));
end
if ~isempty(exinfo.sdfs.extras) % blank
    plot(exinfo.times/10, exinfo.sdfs.extras{1}.sdf, 'r:'); ho
    leg{i+1} = sprintf('blank \t n=%1.0f', exinfo.sdfs.extras{1}.n);
end
ylim_ = get(baseplot, 'ylim');

% window of noise calculation
fill([exinfo.times(201)/10,  exinfo.times(201)/10, ...
    exinfo.times(400)/10, exinfo.times(400)/10] , ...
    [0 ylim_(2)/10 ylim_(2)/10 0], g, 'EdgeColor', g); ho


% annotations
legend(leg, 'Location', 'EastOutside');
s = horzcat(exinfo.sdfs.s{:}); meanfr = mean(mean(s(201:400),2));
title(sprintf('base lat: %1.1f (hmax %1.1f), dur: %1.1f, \n average sd: %1.2f, mean fr: %1.2f',...
    exinfo.lat, exinfo.lat2Hmax, exinfo.dur, mean(sqrt(exinfo.resvars(201:400))), meanfr));
xlim([0 160]); grid on;
ylabel('spk/s');


%------------------------------------------- drug
drugplot = subplot(3,1,2);

% stimulus triggered SDFs 
for i = 1:length(exinfo.sdfs_drug.s) % orientation
    plot(exinfo.times_drug/10, exinfo.sdfs_drug.s{i}, 'Color', c(i,:)); ho
    leg{i} = sprintf('%1.0f \t n=%1.0f', exinfo.sdfs_drug.x(i), ...
        exinfo.sdfs_drug.n(i));
end
if ~isempty(exinfo.sdfs_drug.extras) % blank
    plot(exinfo.times_drug/10, exinfo.sdfs_drug.extras{1}.sdf, 'r:'); ho
    leg{i+1} = sprintf('blank \t n=%1.0f', exinfo.sdfs_drug.extras{1}.n); 
end

% window for noise calculation
fill([exinfo.times_drug(201)/10,  exinfo.times_drug(201)/10, ...
    exinfo.times_drug(400)/10, exinfo.times_drug(400)/10 ], ...
    [0 ylim_(2)/10 ylim_(2)/10 0], g, 'EdgeColor', g); 

legend(leg, 'Location', 'EastOutside');

s_drug = horzcat(exinfo.sdfs_drug.s{:});
meanfr_drug = mean(mean(s_drug(201:400),2));

% annotations
title(sprintf('drug lat: %1.1f (hmax %1.1f), dur: %1.1f, \n average sd: %1.2f, meanfr: %1.2f',...
    exinfo.lat_drug, exinfo.lat2Hmax_drug, exinfo.dur_drug, mean(sqrt(exinfo.resvars_drug(201:400))), ...
    meanfr_drug));
xlim([0 160]); grid on;


%--------------------------------- equalize y axis and plot latency line
ylim_ = [0, max([max(get(baseplot, 'ylim')), max(get(drugplot, 'ylim'))])]; 
set(baseplot, 'ylim', ylim_); % equalize y axis
set(drugplot, 'ylim', ylim_); % equalize y axis

lat = exinfo.lat; lat_drug = exinfo.lat_drug;
lathmax = exinfo.lat2Hmax; lathmax_drug = exinfo.lat2Hmax_drug;
dur = exinfo.dur; dur_drug = exinfo.dur_drug;

plot(baseplot, [lat lat], ylim_, 'k');
plot(baseplot, [lathmax lathmax], ylim_, 'k--');
plot(baseplot, [lat+dur, lat+dur], ylim_, 'k');

plot(drugplot, [lat_drug lat_drug], ylim_, 'k');
plot(drugplot, [lathmax_drug lathmax_drug], ylim_, 'k--');
plot(drugplot, [lat_drug+dur_drug, lat_drug+dur_drug], ylim_, 'k');


%--------------------------------- third plot showing normalized sd(sdfs)

s3 = subplot(3,1,3);
c = getCol(exinfo);

% normalzed deviations from mean response for both conditions
plot(exinfo.times/10, sqrt(exinfo.resvars)./ max(sqrt(exinfo.resvars)), ...
    'Color', lines(1)); hold on;
plot(exinfo.times_drug/10, sqrt(exinfo.resvars_drug)./ ...
    max(sqrt(exinfo.resvars_drug)), c);

plot([lat lat], get(gca, 'YLim'), 'Color', lines(1));
plot([lat_drug lat_drug], get(gca, 'YLim'), c);

% annotations
xlabel('time');
legend('base', exinfo.drugname);
grid on;
xlim([0 160]);

end

function plotRCxORxCO(exinfo)

nplot = length(exinfo.sdfs.y(1,:));

c = hsv(length(exinfo.sdfs.s));
%------------------------------------------ baseline

[~, co_idx]= sort(exinfo.sdfs.y(1,:));

% responses
for co = 1:length(co_idx)
    s(co) = subplot(3, nplot, co);
    
    for i = 1:length(exinfo.sdfs.s)
        plot(exinfo.times/10, exinfo.sdfs.s{i, co_idx(co)}, 'Color', c(i,:), ...
            'DisplayName', sprintf('%1.0f', exinfo.sdfs.n(i,co_idx(co)))); ho
    end
    
    if ~isempty(exinfo.sdfs.extras)
        plot(exinfo.times/10, exinfo.sdfs.extras{1}.sdf, 'r:', ...
            'DisplayName', 'blank'); ho;  
    end
    title(sprintf('CO: %1.2f', exinfo.sdfs.y(1,co_idx(co))));
    l = legend('show','Location', 'EastOutside'); l.FontSize = 6;
    l.Box = 'off';
    
    
     xlim([0 160]);

end

axes('Position', [0.1 0.95 0.8 0.05]);   
off = 0;
for i = 1:length(exinfo.sdfs.s)
    
    plot([0+i+off 1+i+off], [0 0], 'Color', c(i,:)); hold on
    text(0+i+off,0, sprintf('%1.1f', exinfo.sdfs.x(i,1)))
    off = off+0.5;
end
axis off


s(1).YLabel.String = 'Baseline';
set(s, 'YLim', [0 max(max(horzcat(exinfo.sdfs.s{:})))]);
grid on;


%------------------------------------------- drug
[~, co_idx]= sort(exinfo.sdfs_drug.y(1,:));

% responses
for co = 1:length(co_idx)
    s2(co) = subplot(3,nplot,co+nplot);
    
    for i = 1:length(exinfo.sdfs_drug.s)
        plot(exinfo.times_drug/10, exinfo.sdfs_drug.s{i, co_idx(co)}, 'Color', c(i,:), ...
            'DisplayName', sprintf('%1.0f', exinfo.sdfs_drug.n(i,co_idx(co)))); ho
    end
    
    if ~isempty(exinfo.sdfs.extras)
        plot(exinfo.times_drug/10, exinfo.sdfs_drug.extras{1}.sdf, 'r:', ...
            'DisplayName', 'blank'); ho;  
    end
    xlim([0 160]); l = legend('show','Location', 'EastOutside'); l.FontSize = 6; 
        l.Box = 'off';
end
xlabel(exinfo.drugname);
set(s2, 'YLim', [0 max(max(horzcat(exinfo.sdfs_drug.s{:})))]);

%--------------------------------- third plot showing normalized sdf dev


% amplified deviation to mean response
c = getCol(exinfo);

for co = 1:length(co_idx)
    
    s3(co) = subplot(3,nplot,nplot*2+co);
    plot(exinfo.times/10, ...
        sqrt(exinfo.resvars(:,co_idx(co)))./ max(sqrt(exinfo.resvars(:,co_idx(co)))), ...
        'Color', lines(1), 'LineWidth', 1); hold on;
    plot(exinfo.times_drug/10, ...
        sqrt(exinfo.resvars_drug(:,co_idx(co)))./ max(sqrt(exinfo.resvars_drug(:,co_idx(co)))), ...
        c, 'LineWidth', 1);
    title(sprintf('lat: %1.2f (%1.2f) \n lat_drug: %1.2f (%1.2f) ', ...
        exinfo.sdfs.latFP(co_idx(co)), exinfo.sdfs.lat2hmax(co_idx(co)), ...
        exinfo.sdfs_drug.latFP(co_idx(co)), exinfo.sdfs_drug.lat2hmax(co_idx(co))), 'FontSize', 8);
    xlim([0 160]);
end

s2(1).YLabel.String = 'SD(SDF)';

legend('base', exinfo.drugname);
grid on;

 
end