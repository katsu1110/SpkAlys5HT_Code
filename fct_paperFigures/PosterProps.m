function PosterProps(ax, xlim_, ylim_, varargin)
%POSTERPROPS transforms plots to general poster style
%
%
% @author Corinna Lorenz
% @date 8.3.2015
% PosterProps(ax, xlim, ylim)
%
% PosterProps takes the axis and x and y limits and adapts the style. In
% default, this means:
% - scatter data are increased in size to 200 (it does not need to be a
% scatter plot)
% - font size is set to 25
% - transparancy is set to 80%
% - ticks turn outside
% - axis are square
% - box is set off
% - axis tick label are beginning and end of the limit
%
%
%
%
% Using optional arguments, you can change the setting:
% - 'sz' and 'fsz' change scatter data size, respectively font size
%     PosterProps(ax, xlim, ylim, 'sz', 20)
% - 'alpha' is used to change the transperency parameter, must be [0,1]
%     PosterProps(ax, xlim, ylim, 'alpha', 0.5)
% - 'untity' adds a grey unity line and puts it in the background
%     PosterProps(ax, xlim, ylim, 'unity')
% - 'nounity' deletes preexisting, dashed lines counteracting with 'unity'
% - 'save' saves the figure as .svg under the given name
%
% %
% Other setttings that can not be altered:
% - LineWidth is set to 2.
% - Figure position is set. If you need to change it, change it outside,
% after you call PosterProps
%
%
% NOTE: This function was written under use of Matlab R2015a. There have
% been major changes in the plotting functions between 2014 to 2015,
% especially with the transparancy. I would recommend to use R2015a or
% later versions to work with PosterProps (also because I think plots look
% nicer with new versions, otherwise try to make a copy and adapt it to
% older versions.
%
% NOTE2: Put the unity function in the same folder to guarantee smooth
% collaborations.

h = setPaperPlotsProps();delete(h);
global uni_k uni_r;


% set(gcf, 'Position', [492   269   550   500]);
box off

fname = '';

fsz = 25;
sz = 200;
transpy = 0.5;
unity_flag = false;
cross_flag = false;
square_flag = false;
j=1;
while j <= length(varargin)
    switch varargin{j}
        case 'sz'
            sz = varargin{j+1};
            j=j+1;
        case 'fsz'
            fsz = varargin{j+1};
            j=j+1;
        case 'unity'
            unity_flag = true;
        case 'cross'
            cross_flag = true;
        case 'alpha'
            transpy = varargin{j+1};
            j=j+1;
        case 'square'
            square_flag= true;
        case 'save'
            fname = varargin{j+1};
            j=j+1;
        otherwise
            j=j+1;
    end
    j = j+1;
end


set(ax, 'xlim', [xlim_(1) xlim_(end)], ...
    'ylim', [ylim_(1) ylim_(end)], ...
    'XTick', xlim_, 'YTick', ylim_, ...
    'TickDir', 'out', ...
    'LineWidth', 0.5, ...
    'XTickLabel', cellfun(@num2str, num2cell(xlim_'),'UniformOutput', 0), ...
    'YTickLabel', cellfun(@num2str, num2cell(ylim_'),'UniformOutput', 0), ...
    'LineWidth', 0.5, ...
    'FontSize', fsz, ...
    'XColor', 'k', 'YColor', 'k');
%     'Position', [0.2 0.2 0.6 0.6],


% format scatter plot
set(findobj(ax, 'Type', 'Scatter'), 'SizeData', sz, 'MarkerFaceAlpha', transpy, ...
 'MarkerEdgeAlpha', 0.6, 'MarkerEdgeColor', 'w', 'LineWidth', 0.3); 


% format line plot with specified data
set(findobj(ax, 'Type', 'Line', 'LineStyle', 'o'), 'MarkerSize',sz, ...
    'MarkerFaceAlpha', transpy, 'MarkerEdgeAlpha', 0.3);

% format error bar
set(findobj(ax,'Type', 'Errorbar'), 'MarkerSize', 2, 'LineWidth', 0.35);

% format dashed lines
set(findobj(ax,'Type', 'Line', 'LineStyle', '--'), 'LineStyle', ':');


% change black color to dark grey
set(findobj(ax, 'Type', 'Line', 'Color', 'k'),'Color',uni_k);
set(findobj(ax, 'Type', 'Scatter', 'MarkerFaceColor', 'k'),'MarkerFaceColor',uni_k);
set(findobj(ax, 'Type', 'Errorbar', 'Color', 'k'), 'Color',uni_k);

set(findobj(ax, 'Type', 'Line', 'Color', 'r'),'Color',uni_r);
set(findobj(ax, 'Type', 'Scatter', 'MarkerFaceColor', 'r'),'MarkerFaceColor',uni_r);
set(findobj(ax, 'Type', 'Errorbar', 'Color', 'r'), 'Color',uni_r);


% add unity or cross. note that both are not compatible with each other 
if unity_flag; delete(findobj(ax, 'type','Line')); unity; end
if cross_flag; delete(findobj(ax, 'type','Line')); crossl;end
if square_flag; axis square; end


ax.XLabel.String = ''; 
ax.YLabel.String = '';
ax.Layer = 'bottom';

% print to svg file if wanted
if ~isempty( fname )
    print(gcf, fname, '-dsvg');
end

