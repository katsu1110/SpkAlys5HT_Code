function plotTC_compRC( info_ )



tit = sprintf(['#%1.0f ' info_.cluster  ...
    ' me=%1.0f, RC=true \n'  ...
    ' gain=%1.2f add= %1.2f \n'...
    ' r2drug=%1.3f    r2base=%1.3f   r2both=%1.3f '],...
    info_.id, ...
    info_.ocul, info_.gslope, info_.yoff, ...
    info_.rsqr_drug, info_.rsqr_cont, info_.rsqr_both);

%% tuning curve plot
parvals = info_.ratepar;
parvals_drug = info_.ratepar_drug;


subplot(1,2,1)
pos = get(gca, 'OuterPosition');
set(gca, 'OuterPosition', [pos(1:3), 0.9])


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
    info_.ratemn, ...
    info_.ratesme, ...
    'Color', c); hold on;

errorbar(1:length(parvals_drug), ...
    info_.ratemn_drug, ...
    info_.ratesme_drug, ...
    'LineStyle', ':' , 'Color', c);

legend('base', info_.drugname);

parvalab = cellstr(num2str(parvals')); 
set(gca, 'XTick', 1:length(parvals), 'XTickLabel', parvalab);
xlim([0, length(parvalab)+1]); xlabel(info_.param1), ylabel('mean spike rate (\pm std)');



%% gain change plot
subplot(1,2,2);
pos = get(gca, 'OuterPosition');
set(gca, 'OuterPosition', [pos(1:3), 0.9]);

cont = info_.ratemn;
drug = info_.ratemn_drug;

scatter(cont, drug, 'MarkerEdgeColor', 'k'); hold on;
plot(cont, cont.*info_.gslope + info_.yoff, 'k'); hold on;
    
    
xlim([ min([cont, drug]), max([cont, drug])+0.001 ]);
ylim([ min([cont, drug]), max([cont, drug])+0.001 ]);

xlabel('control [mean spk/s]');
ylabel([info_.drugname '[mean spk/s]']);

%% title


ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible', 'off','Units','normalized', 'clipping', 'off');
text(0.5, 1, tit,'HorizontalAlignment' ,'center','VerticalAlignment', 'top')


end

