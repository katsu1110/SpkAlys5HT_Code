function plotTC_comp( Trials, Trials_drug, info_ )



tit = sprintf(['#%1.0f ' info_.cluster  ...
    ' me=%1.0f'  ...
    ' RC=%1.0f'  ' adapt=%1.0f' '\n ' ...
    ' gain=%1.2f add= %1.2f \n'...
    ' r2drug=%1.3f    r2base=%1.3f   r2both=%1.3f '],...
    info_.id, ...
    info_.ocul, info_.isRC, info_.isadapt,info_.gslope,  info_.yoff, ...
    info_.rsqr_drug, info_.rsqr_cont, info_.rsqr_both);

%% tuning curve plot
parvals = unique([Trials.(info_.param1)]);

subplot(1,2,1)
pos = get(gca, 'OuterPosition');
set(gca, 'OuterPosition', [pos(1:3), 0.9])

for i = 1:length(parvals)
    errbar{i}       = [Trials(  parvals(i) == [Trials.(info_.param1)] ).spkRate];    
    errbar_drug{i}  = [Trials_drug(  parvals(i) == [Trials_drug.(info_.param1)]  ).spkRate];
end

switch info_.ocul
    case -1
        c = [0 0 1];
    case 0
        c = [0 0.5 0];
    case 1
        c = [1 0 0];
end

% plot errorbars
errorbar(1:length(parvals), ...
    cellfun(@mean, errbar), ...
    cellfun(@std, errbar)/sqrt(length(errbar)),...
    'Color', c); hold on;

errorbar(1:length(parvals), ...
    cellfun(@mean, errbar_drug), ...
    cellfun(@std, errbar_drug)/sqrt(length(errbar_drug)), ...
    'LineStyle', ':' , 'Color', c);


legend('base', info_.drugname);

% plot number of data points
for kk = 1:length(parvals)
    text(kk, mean(errbar{kk}),...
        num2str(length(errbar{kk})), 'Color', c);
    
    text(kk, mean(errbar_drug{kk}),...
        num2str(length(errbar_drug{kk})), 'Color', c);
end
hold off;


parvalab = cellstr(num2str(parvals')); 
parvalab{end} = 'blank';
set(gca, 'XTick', 1:length(parvals), 'XTickLabel', parvalab);
xlim([0, length(parvalab)+1]); xlabel(info_.param1), ylabel('mean spike rate (\pm std)');

axis square;
box off;


%% gain change plot
subplot(1,2,2);
pos = get(gca, 'OuterPosition');
set(gca, 'OuterPosition', [pos(1:3), 0.9]);


cont = info_.ratemn (ismember(info_.ratepar, info_.ratepar_drug)) ;
drug = info_.ratemn_drug (ismember(info_.ratepar_drug, info_.ratepar)) ;


scatter(cont, drug, 'MarkerEdgeColor', 'k'); hold on;

if info_.is5HT
plot([min(get(gca, 'xlim')) max(get(gca, 'xlim'))] ,...
    [ min(get(gca, 'xlim'))*info_.gslope + info_.yoff, ...
    max(get(gca, 'xlim')).*info_.gslope + info_.yoff], 'r'); hold on;
else
    plot([min(get(gca, 'xlim')) max(get(gca, 'xlim'))] ,...
    [ min(get(gca, 'xlim'))*info_.gslope + info_.yoff, ...
    max(get(gca, 'xlim')).*info_.gslope + info_.yoff], 'k'); hold on;
end
unity
    
xlabel('control [mean spk/s]');
ylabel([info_.drugname '[mean spk/s]']);

axis square;
box off;

%% title


ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible', 'off','Units','normalized', 'clipping', 'off');
text(0.5, 1, tit,'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
axis square;
box off;


end

