function setHist(ax, lim_, varargin)
% set the histogram for paper plotting


global uni_k uni_r;

ylim_h = 50;
arrowpos = [0 0];
fsz = 8; msz = 2; txtsz = 5; 
transpy = 0.8;
j=1; rot = 0;
while j <= length(varargin)
    switch varargin{j}
        case 'ylim_h'
            ylim_h = varargin{j+1};
        case 'arrow'
            arrowpos = varargin{j+1};
        case 'fsz'
            fsz = varargin{j+1};
        case 'txtsz'
            txtsz = varargin{j+1};
        case 'sz'
            sz = varargin{j+1};
        case 'alpha'
%             transpy = varargin{j+1};
        case 'msz'
            msz = varargin{j+1};
        case 'rot'
            rot = varargin{j+1};
    end
    j = j+2;
end

PosterProps(ax, lim_, [0 ylim_h], varargin{:}); 

ax.YTick = ylim_h;
ax.YTickLabel = ax.YTickLabel(2);
ax.XTick = [];
ax.FontSize = fsz;
delete(ax.Children(1));
delete(ax.Children(2));
ax.Children(1).XData = arrowpos(1);
ax.Children(1).YData = ax.YTick+(ax.YTick*0.1) ;
ax.Children(1).Color = uni_r;
ax.Children(1).MarkerFaceColor = uni_r;
ax.Children(1).MarkerSize = msz;

ax.Children(2).XData = arrowpos(2);
ax.Children(2).YData = ax.YTick+(ax.YTick*0.1) ;
ax.Children(2).Color = uni_k;
ax.Children(2).MarkerFaceColor = uni_k;
ax.Children(2).MarkerSize = msz;
ax.Title.String = '';

if ax.View(1) == 90
    ax.Children(1).Marker = '<';
    ax.Children(2).Marker = '<';
end

% adjust bin edges and width
if ~strcmp(ax.Children(3).BinWidth, 'nonuniform')
    if ax.Children(3).BinWidth < ax.Children(4).BinWidth
        ax.Children(4).BinWidth = ax.Children(3).BinWidth;
    else
        ax.Children(3).BinWidth = ax.Children(4).BinWidth;
    end
end

% set transparancy
ax.Children(3).FaceAlpha = transpy; ax.Children(4).FaceAlpha = transpy;
ax.Children(3).FaceColor = uni_k; ax.Children(4).FaceColor = uni_r;
ax.Children(3).EdgeAlpha = 0.2; ax.Children(4).EdgeAlpha = 0.2;
ax.Children(3).EdgeColor = 'w'; ax.Children(4).EdgeColor = 'w';

% add statistic mean/geometric mean/median as text
text(ax, arrowpos(1)-(ax.XLim(2)*0.01), ax.YTick+(ax.YTick*0.4), num2str(arrowpos(1)), ...
    'Color', uni_r, 'FontSize', txtsz, 'Rotation', rot, 'HorizontalAlignment', 'center');
text(ax, arrowpos(2)+(ax.XLim(2)*0.01), ax.YTick+(ax.YTick*0.4), num2str(arrowpos(2)),...
    'Color', uni_k, 'FontSize', txtsz, 'Rotation', rot, 'HorizontalAlignment', 'center');

ax.Clipping = 'off';
% axis off;
end
