function PosterProps_wHist(fig, xlim_, ylim_, varargin)
%POSTERPROPS transforms plots to general poster style
%
%
% @author Corinna Lorenz
% @date 06.12.2016
% PosterProps_wHist(ax, xlim_, ylim_, varargin)
%
% PosterProps takes the axis and x and y limits and adapts the axis
% accordingly also for the histograms.
%
%
% Using optional arguments, you can adjust the following:
% - 'ylim_uph' + argument y-limit of the upper histogram
% - 'ylim_rih' same for the right one
% - 'arrow_uph' two element vector containing the x-position of arrows in
% the upper histogram  
% - 'arrow_rih' same for the right histogramm
% 
% all other arguments are given to PosterProps.
% 
% 
% 
% NOTE: This function was written under use of Matlab R2016a. There have
% been major changes in the plotting functions between 2014 to 2015,
% especially with the transparancy. I would recommend to use R2015a or
% later versions to work with PosterProps (also because I think plots look
% nicer with new versions, otherwise try to make a copy and adapt it to
% older versions.


fname = '';
ylim_uph = 50;
ylim_rih = 50;
arrow_uph = [0 0];
arrow_rih = [0 0];
j=1;
while j <= length(varargin)
    switch varargin{j}
        case 'ylim_uph'
            ylim_uph = varargin{j+1};
        case 'ylim_rih'
            ylim_rih = varargin{j+1};
        case 'arrow_uph'
            arrow_uph = varargin{j+1};
        case 'arrow_rih'
            arrow_rih = varargin{j+1};
        case 'saveeps'
            fname = varargin{j+1};
    end
    j = j+2;
end


ax = get(fig, 'Children');

axmain = ax(1);
axhist1 = ax(2);
axhist2 = ax(3);

delete(findobj(axmain, 'LineStyle', '-'));

%%% main plot
PosterProps(axmain, xlim_, ylim_, varargin{:}, 'cross'); 
set(axmain, 'Position', [0.2 0.2 0.5 0.5])
axmain.Title.String = '';


%%% right histogram
PosterProps(axhist1, ylim_, [0 ylim_rih], varargin{:}); 

axhist1.YTick = ylim_rih;
axhist1.YTickLabel = axhist1.YTickLabel(2);
axhist1.XTick = [];
axhist1.Position = [0.7 0.2 0.15 0.5];
delete(axhist1.Children(1));
delete(axhist1.Children(2));
axhist1.Children(1).XData = arrow_rih(1);
axhist1.Children(1).YData = axhist1.YTick;
axhist1.Children(1).Color = 'r';
axhist1.Children(1).MarkerFaceColor = 'r';
axhist1.Children(1).Marker = '<';

axhist1.Children(2).XData = arrow_rih(2);
axhist1.Children(2).YData = axhist1.YTick;
axhist1.Children(2).Marker = '<';
axhist1.Title.String = '';

if axhist1.Children(3).BinWidth < axhist1.Children(4).BinWidth
    axhist1.Children(4).BinEdges = axhist1.Children(3).BinEdges;
else
    axhist1.Children(3).BinEdges = axhist1.Children(4).BinEdges;
end

text(axhist1, arrow_rih(1)-(axhist1.XLim(2)*0.1), axhist1.YTick+(axhist1.YTick*0.1), num2str(arrow_rih(1)));
text(axhist1, arrow_rih(2)+(axhist1.XLim(2)*0.1), axhist1.YTick+(axhist1.YTick*0.1), num2str(arrow_rih(2)));


%%% upper histogram
PosterProps(axhist2, xlim_, [0 ylim_uph],varargin{:}); 
axhist2.YTick = ylim_uph;
axhist2.YTickLabel = axhist2.YTickLabel(2);
axhist2.XTick = [];
axhist2.Position = [0.24 0.73 0.422 0.15];
delete(axhist2.Children(1));
delete(axhist2.Children(2));
axhist2.Children(1).XData = arrow_uph(1);
axhist2.Children(1).YData = axhist2.YTick;
axhist2.Children(1).Color = 'r';
axhist2.Children(1).MarkerFaceColor = 'r';

axhist2.Children(2).XData = arrow_uph(2);
axhist2.Children(2).YData = axhist2.YTick;
axhist2.Title.String = '';

if axhist2.Children(3).BinWidth < axhist2.Children(4).BinWidth
    axhist2.Children(4).BinEdges = axhist2.Children(3).BinEdges;
else
    axhist2.Children(3).BinEdges = axhist2.Children(4).BinEdges;
end

text(axhist2, arrow_uph(1)-(axhist2.XLim(2)*0.02), axhist2.YTick+(axhist2.YTick*0.1), num2str(arrow_uph(1)));
text(axhist2, arrow_uph(2)+(axhist2.XLim(2)*0.02), axhist2.YTick+(axhist2.YTick*0.1), num2str(arrow_uph(2)));


set(fig, 'Position', [958   380   708   598]);


if ~isempty( fname )
    print(gcf, fname, '-dpdf');
end



function ticks = getTicks(scl, lim)

if strcmp(scl, 'log')
    ticks = [lim(1) 1 lim(2)] ;
else
    ticks = [lim(1) 0 lim(2)] ;
end

 