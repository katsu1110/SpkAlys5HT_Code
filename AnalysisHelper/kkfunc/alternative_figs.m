function alternative_figs(nc, varargin)

ms = 40;
lw = 0.5;

if nargin==0
    nc = 0;
end

close all

if nc==1
    % noise correlation
    uiopen('C:\Users\katsuhisa\Documents\code\analysis\integrated\interaction_project\SpkAlys5HT_Code\AnalysisHelper\kkfunc\nc_8.fig',1)
    label = 'Noise Correlation';
    axrange = [-1 1];
    axlab = [-1 0 1];
elseif nc==0
    uiopen('C:\Users\katsuhisa\Documents\code\analysis\integrated\interaction_project\SpkAlys5HT_Code\AnalysisHelper\kkfunc\ff_8.fig',1)
    label = 'Fano Factor';
    axrange = [0 11];
    axlab = [0 10];
end

% current figure handle;
h = gcf;

% axis object
axesObjs = get(h, 'Children'); 

% data object inside the axis
dataObjs = get(axesObjs, 'Children'); 

% scatter and line
sdata = findobj(dataObjs, 'Type', 'scatter');
ldata = findobj(dataObjs, 'Type', 'line');

% color code
col_5HT = hot(18);
col_NaCl = gray(12);

% split data based on stimuli
animal = ones(1, length(sdata));
drug = zeros(1, length(sdata));
sig = ones(1, length(sdata));
stmdim = zeros(1, length(sdata));
for i = 1:length(sdata)
    % animal
    if strcmp(sdata(i).Marker, 'square')
        animal(i) = 0;
    end
    
    % significance
    if isequal(sdata(i).MarkerFaceColor, ones(1,3))
        sig(i) = 0;
    end
    
    % drug & stimulus
    if isequal(sdata(i).MarkerEdgeColor, col_5HT(4,:))
        drug(i) = 1;
    elseif isequal(sdata(i).MarkerEdgeColor, col_5HT(6,:))
        drug(i) = 1;
        stmdim(i) = 1;
    elseif isequal(sdata(i).MarkerEdgeColor, col_NaCl(3,:))
        stmdim(i) = 1;    
    elseif isequal(sdata(i).MarkerEdgeColor, col_5HT(8,:))
        drug(i) = 1;
        stmdim(i) = 2;
    elseif isequal(sdata(i).MarkerEdgeColor, col_NaCl(5,:))
        stmdim(i) = 2;
    elseif isequal(sdata(i).MarkerEdgeColor, col_5HT(10,:))
        drug(i) = 1;
        stmdim(i) = 3;
    elseif isequal(sdata(i).MarkerEdgeColor, col_NaCl(7,:))
        stmdim(i) = 3;
    end
end

% replot OR data
figure;
for i = 1:length(ldata)
    plot(ldata(i).XData, ldata(i).YData, '-', 'color',...
        ldata(i).Color, 'lineWidth', ldata(i).LineWidth)
    hold on;
end

for i = 1:length(sdata)-1
    if stmdim(i)==0
        % edge color
        if drug(i)==1
            ecol = 'r';
        elseif drug(i)==0
            ecol = 'k';
        end

        % face color
        if sig(i)==0
            fcol = ones(1,3);
        elseif sig(i)==1
            fcol = ecol;
        end

        % scatter
        hold on;
        scatter(sdata(i).XData, sdata(i).YData, ms,...
            sdata(i).Marker, 'linewidth', lw,...
            'markerfacecolor', fcol, 'markeredgecolor', ecol,...
            'markerfacealpha', sdata(i).MarkerFaceAlpha,...
            'markeredgealpha', sdata(i).MarkerEdgeAlpha)
    end
end

% cosmetics
title(label)
xlabel('baseline')
ylabel('5HT / NaCl')
xlim(axrange)
ylim(axrange)
set(gca, 'XTick', axlab)
set(gca, 'YTick', axlab)
set(gca,'box','off')
set(gca,'TickDir','out')
set(gca,'FontSize',16)
axis square

% bar graph
barmat = nan(2,4);

for i = 1:2
    basedata = [sdata(drug==i-1).XData];
    drugdata = [sdata(drug==i-1).YData];
    
%     basedata = [sdata(drug==i-1 & stmdim==0).XData];
%     drugdata = [sdata(drug==i-1 & stmdim==0).YData];
    
    nans = isnan(basedata) | isnan(drugdata);
    basedata(nans) = [];
    drugdata(nans) = [];
    
    % baseline mean
    barmat(i, 1) = mean(basedata);
%     barmat(i, 1) = median(basedata);
    
    % drug mean
    barmat(i,2) = mean(drugdata);
%     barmat(i,2) = median(drugdata);
    
    % baseline SEM
    barmat(i,3) = std(basedata)/sqrt(length(basedata));
    
    % drug SEM
    barmat(i,4) = std(drugdata)/sqrt(length(drugdata));
end

figure;
b = bar(1, barmat(1,1),0.975);
b.FaceColor = 'w';
b.EdgeColor = 'r';
b.FaceAlpha = 0.5;
b.EdgeAlpha = 0.8;
hold on;
b = bar(2, barmat(1,2),0.975,'r');
b.EdgeColor = 'w';
b.FaceAlpha = 0.5;
hold on;
b = bar(4, barmat(2,1),0.975,'w');
b.EdgeColor = 'k';
b.FaceAlpha = 0.5;
b.EdgeAlpha = 0.8;
hold on;
b = bar(5, barmat(2,2),0.975,'k');
b.EdgeColor = 'w';
b.FaceAlpha = 0.5;
hold on;
e = errorbar(1:2, barmat(1,1:2), barmat(1,3:4));
e.LineStyle = 'none';
e.MarkerSize = 10;
e.Color = 'r';
e.CapSize = 0;
hold on;
e = errorbar(4:5, barmat(2,1:2), barmat(2,3:4));
e.LineStyle = 'none';
e.MarkerSize = 10;
e.Color = 'k';
e.CapSize = 0;
hold on;
plot([0 6],[0 0],'-','color',zeros(1,3),'linewidth',0.25)

% % individual data
% for i = 1:length(sdata)
%         % edge color & range
%         if drug(i)==1
%             ecol = 'r';
%             xx = [1 2];
%         elseif drug(i)==0
%             ecol = 'k';
%             xx = [4 5];
%         end
% 
%         % face color
%         if sig(i)==0
%             fcol = ones(1,3);
%         elseif sig(i)==1
%             fcol = ecol;
%         end
%             
%         hold on;
%     
%         for k = 1:length(sdata(i).XData)
%             scatter(xx, [sdata(i).XData(k), sdata(i).YData(k)], ms,...
%                     sdata(i).Marker, 'linewidth', lw,...
%                     'markerfacecolor', fcol, 'markeredgecolor', ecol,...
%                     'markerfacealpha', sdata(i).MarkerFaceAlpha,...
%                     'markeredgealpha', sdata(i).MarkerEdgeAlpha)
%                 hold on;
%         end
% end


% cosmetics
yy = get(gca, 'YLim');
% set(gca, 'XTick', [1.5 4.5], 'XTickLabel', {'5HT', 'NaCl'})
set(gca, 'XTick', [1 2 4 5], 'XTickLabel', {'baseline', '5HT', 'baseline', 'NaCl'})
xtickangle(45)
set(gca, 'YTick', [0 yy(2)]);
ylabel(label)
set(gca,'box','off')
set(gca,'TickDir','out')
set(gca,'FontSize',16)
axis square



