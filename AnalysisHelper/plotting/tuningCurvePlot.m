function h = tuningCurvePlot(exinfo)
% h = tuningCurvePlot(exinfo)
% 
% h is the figure handle of the figure containing the superimposed raw and
% fitted tuning curves for both experiments +/- sme using the information
% from exinfo.fitparam and exinfo.fitparam_drug. the data is provided by
% evalSingleDG where, depending on the stimulus feature, descriptive
% functions were fitted to the raw tuning curves.
% 
% the figure is saved as exinfo.fig_tc
% 
% @CL

h = figure('Name', exinfo.figname);


if strcmp(exinfo.param1, 'or') 
    fittedTC_or(exinfo)
elseif  strcmp(exinfo.param1, 'sf')
    
    s1 =subplot(2,1,1); % linearly scaled and fitted tc
    fittedTC_sf(exinfo, exinfo.fitparam.others{1}, exinfo.fitparam_drug.others{1});
    s2 = subplot(2,1,2); % log scaled and fitted tc
    fittedTC_sf(exinfo, exinfo.fitparam.others{2}, exinfo.fitparam_drug.others{2});    
    
elseif strcmp(exinfo.param1, 'co')
    fittedTC_co(exinfo);
elseif strcmp(exinfo.param1, 'sz')
    fittedTC_sz(exinfo);
end

set(h, 'UserData', exinfo);
savefig(h, exinfo.fig_tc);
close(h);

end



%% size data
function fittedTC_sz(exinfo)
% size tuning curve --> see fitSZ()


c = getCol(exinfo);
val = exinfo.fitparam.val;
errorbar( val.sz,val.mn, val.sem, 'o', 'Color', c, 'MarkerFaceColor', c); ho
text( val.sz,val.mn, num2str(exinfo.nrep(exinfo.ratepar<1000)));
plot( val.x, val.y, 'Color', c);

val = exinfo.fitparam_drug.val;
errorbar( val.sz,val.mn, val.sem, 'o', 'Color', c, 'MarkerFaceColor', 'w');
plot( val.x, val.y, 'Color', c, 'LineStyle', '--');
text( val.sz,val.mn, num2str(exinfo.nrep_drug(exinfo.ratepar_drug<1000)));

plotSpontResp(exinfo, c);
% plot( ones(2,1)*exinfo.fitparam.mu, get(gca, 'YLim'), 'k');
% plot( ones(2,1)*exinfo.fitparam_drug.mu, get(gca, 'YLim'), 'k--');

ylim_ = get(gca, 'YLim');
ylim([0 ylim_(2)]);

title(sprintf(['Base: pref sz : %1.3f, SI= %1.3f, r2=%1.3f\n' exinfo.drugname ...
    'pref sz: %1.3f, SI=%1.3f, r2=%1.3f'],...
    exinfo.fitparam.mu, exinfo.fitparam.SI, exinfo.fitparam.r2, ...
    exinfo.fitparam_drug.mu, exinfo.fitparam_drug.SI, exinfo.fitparam_drug.r2));
xlabel('size'), ylabel('spks/s');

end


%% orientation contrast RC
function fittedTC_orco(exinfo)
% 2D tuning of contrast and orientation 

n = ceil(length(exinfo.fitparam.others.OR)/2)+1;
for i = 1:length(exinfo.fitparam.others.OR)
    subplot( n , 2, i)
    exinfo_temp = exinfo;
    exinfo_temp.fitparam  = exinfo.fitparam.others.OR(i);    
    exinfo_temp.fitparam_drug  = exinfo.fitparam_drug.others.OR(i);
    fittedTC_or(exinfo_temp); xlabel(''); ylabel('');
    ylabel(['co = ' num2str(exinfo.sdfs.y(1,i))]);
end


% difference
s_diff = subplot( n, 2, 6);
[~, idx] = sort(exinfo.sdfs.y(1,:));
surf(exinfo.sdfs.x(:,1), exinfo.sdfs.y(1,idx), ...
    exinfo.sdfs.mn_rate(1).mn(:, idx)'- exinfo.sdfs_drug.mn_rate(1).mn(:, idx)'); ho
title('Difference (Base-Drug)');ylabel('co'); xlabel('or'); zlabel('spk/frame');
colormap jet;
set(s_diff, 'YScale', 'log', ...
    'XLim', [0 180], ...
    'YLim', [min(exinfo.sdfs_drug.y(1,:)) max(exinfo.sdfs_drug.y(1,:))])
colorbar



% color axis specification
minc = min(min([exinfo.sdfs.mn_rate(1).mn; exinfo.sdfs_drug.mn_rate(1).mn]));
maxc = max(max([exinfo.sdfs.mn_rate(1).mn; exinfo.sdfs_drug.mn_rate(1).mn]));


% Baseline
s(1) = subplot( n, 2, 7);
[~, idx] = sort(exinfo.sdfs.y(1,:));
surf(exinfo.sdfs.x(:,1), exinfo.sdfs.y(1,idx), exinfo.sdfs.mn_rate(1).mn(:, idx)'); ho
title('Baseline');ylabel('co'); xlabel('or'); zlabel('spk/frame');
caxis([minc maxc]);


% 5HT/NaCl
s(2) = subplot( n, 2, 8);
[~, idx] = sort(exinfo.sdfs_drug.y(1,:));
surf(exinfo.sdfs_drug.x(:,1), exinfo.sdfs_drug.y(1,idx),exinfo.sdfs_drug.mn_rate(1).mn(:, idx)');
title(exinfo.drugname); ylabel('co'); xlabel('or'); zlabel('spk/frame');
caxis([minc maxc]);


set(s, 'YScale', 'log', 'ZLim', ...
    [min([s(1).ZLim(1), s(2).ZLim(1)]),  max([s(1).ZLim(2), s(2).ZLim(2)])], ...
    'XLim', [0 180], 'YLim', [min(exinfo.sdfs_drug.y(1,:)) max(exinfo.sdfs_drug.y(1,:))]);


set(findobj('Type', 'Axes'), 'FontSize', 6)
set(gcf, 'Position', [1229         206         436         763]);
end

%% orientation data
function fittedTC_or(exinfo)
% orientation tuning curve  --> see fitOR()

c = getCol(exinfo); % experiment color
ft = exinfo.fitparam; % fit parameters for the baseline experiment
ft_drug = exinfo.fitparam_drug;% fit parameters for the drug experiment

args_base = {'o', 'Color', c, 'MarkerSize', 5, 'MarkerFaceColor', c};
args_drug = args_base; args_drug{end} = 'w';

fittedTC_or_Helper(ft, args_base, {c}); ho;
fittedTC_or_Helper(ft_drug,  args_drug, {c, 'LineStyle', '--'});

plotSpontResp(exinfo, c)% plot spontaneous response as dotted line

% axes labels
xlabel('orientation'); 
if exinfo.isRC
    ylabel('spks/frame (\pm sme)');
else
    ylabel('spks/s (\pm sme)');
end



% axis settings
minmu = min([ft.mu, ft_drug.mu]);
maxmu = max([ft.mu, ft_drug.mu]);
xlim([minmu-110 maxmu+110]);

title( sprintf( ['B: pf=%1.0f, bw=%1.0f, amp=%1.1f, off=%1.1f, r2=%1.1f \n' ...
    exinfo.drugname ': pf=%1.0f, bw=%1.0f, amp=%1.1f, off=%1.1f r2=%1.1f \n'], ...
    ft.mu, ft.sig, ft.a, ft.b, ft.r2, ...
    ft_drug.mu, ft_drug.sig, ft_drug.a, ft_drug.b, ft_drug.r2), 'FontSize', 8);
end

function fittedTC_or_Helper(ft, errArgs, lineArgs)
    errorbar(ft.val.or, ft.val.mn, ft.val.sem, errArgs{:}); ho
    plot(ft.x, ft.y, lineArgs{:}); ho
end

%% spatial frequency data
function fittedTC_sf(exinfo, ft, ft_drug)
% spatial frequency tuning curve --> see fitSF()


c = getCol(exinfo);
args_base = {'o', 'Color', c, 'MarkerSize', 5, 'MarkerFaceColor', c};
args_drug = args_base; args_drug{end} = 'w';

fittedTC_sf_Helper(ft,  args_base, {c}, exinfo.nrep(exinfo.ratepar<1000)); ho;
fittedTC_sf_Helper(ft_drug,  args_drug, {c, 'LineStyle', '--'}, exinfo.nrep_drug(exinfo.ratepar_drug<1000));

plotSpontResp(exinfo, c)
xlabel('spatial frequency');    ylabel('spks/s (\pm sme)');
axis square; box off;

xdata = [ft.val.sf, ft_drug.val.sf];
if any(xdata<0)
    xlim([min(xdata)-1, max(xdata)+2]);
else
    xlim([0, max(xdata)+2]);
end


title( sprintf( ['B: pf=%1.2f, bw=%1.2f, amp=%1.1f, off=%1.1f r2=%1.3f\n' ...
    exinfo.drugname ': pf=%1.2f, bw=%1.2f, amp=%1.1f, off=%1.1f r2=%1.3f\n'], ...
    ft.mu, ft.sig, ft.a, ft.b, ft.r2, ...
    ft_drug.mu, ft_drug.sig, ft_drug.a, ft_drug.b, ft_drug.r2), 'FontSize', 8);

end


function fittedTC_sf_Helper(ft, errArgs, lineArgs, nrep)

errorbar(ft.val.sf, ft.val.mn, ft.val.sem, errArgs{:}); ho
plot(ft.x, ft.y, lineArgs{:}); ho
text(ft.val.sf, ft.val.mn, num2str(nrep), 'FontSize', 8);

set(gca, 'XLim', [min(ft.val.sf) max(ft.val.sf)]);
end

%% contrast data
function fittedTC_co(exinfo)
% contrast tuning curve --> see fitCO()

c = getCol(exinfo);

mn0 = exinfo.ratemn;            mn1 = exinfo.ratemn_drug;
sme0 = exinfo.ratesme;          sme1 = exinfo.ratesme_drug;   
par0 = exinfo.ratepar;          par1 = exinfo.ratepar_drug;
nrep0 = exinfo.nrep;           nrep1 = exinfo.nrep_drug;
errArgs = {'o', 'Color', c, 'MarkerSize', 5};


idx = par0==0 | par0>1 ;
if any(idx); plot([eps 100], [mn0(idx) mn0(idx)], c); hold on; end;
idx1 = par1==0 | par1>1;
if any(idx); plot([eps 100], [mn1(idx1) mn1(idx1)], c, 'LineStyle', '--'); end;

% plot errorbars
fittedTC_co_Helper(exinfo.fitparam, {c})
fittedTC_co_Helper(exinfo.fitparam_drug, {c, 'LineStyle', '--'})

errorbar(par0.*100, mn0, sme0, errArgs{:}, 'MarkerFaceColor', c); hold on;
text(par0.*100, mn0, num2str(nrep0));
errorbar(par1.*100, mn1, sme1, errArgs{:}, 'MarkerFaceColor', 'w');
text(par1.*100, mn1, num2str(nrep1));

legend('base', exinfo.drugname, 'Location', 'northwest');
set(gca, 'XScale','log', 'XLim', [1 100]);
box off;



title( sprintf(['r_{max} * (c^n / (c^n + c50^n) ) + m \n '...
    'baseline: r_{max}=%1.0f, n=%1.1f, c50=%1.1f, m=%1.1f r2=%1.2f \n' ...
    exinfo.drugname ' ind : r_{max}=%1.0f, n=%1.1f, c50=%1.1f, m=%1.1f r2=%1.2f \n' ...
    '   ' exinfo.drugname ' basepar: r_{max}=%1.0f, n=%1.1f, c50=%1.1f, m=%1.1f r2=%1.2f  \n' ....
    'a_{cog} = %1.1f (%1.2f)  a_{actg} = %1.1f (%1.2f)  a_{resg} = %1.1f (%1.2f) '],...
    exinfo.fitparam.rmax, ...
    exinfo.fitparam.n, exinfo.fitparam.c50, ...
    exinfo.fitparam.m, exinfo.fitparam.r2, ...
    exinfo.fitparam_drug.rmax, ...
    exinfo.fitparam_drug.n, exinfo.fitparam_drug.c50, ...
    exinfo.fitparam_drug.m, exinfo.fitparam_drug.r2, ...
    exinfo.fitparam_drug.rmax, ...
    exinfo.fitparam_drug.sub.n, exinfo.fitparam_drug.sub.c50, ...
    exinfo.fitparam_drug.sub.m, exinfo.fitparam_drug.sub.r2, ...
    exinfo.fitparam_drug.sub.a_cg, exinfo.fitparam_drug.sub.r2_cg, ...
    exinfo.fitparam_drug.sub.a_ag, exinfo.fitparam_drug.sub.r2_ag, ...
    exinfo.fitparam_drug.sub.a_rg, exinfo.fitparam_drug.sub.r2_rg), ...
    'FontSize', 7);
xlabel('contrast (log)');

end

function fittedTC_co_Helper(ft, lineArgs)

plot(ft.x .*100, ft.y, lineArgs{:}); ho

end


%% helper

function plotSpontResp(exinfo, c)
% plot the response to blanks as horizontal line

idx0 = exinfo.ratepar>1000;
idx2 = exinfo.ratepar_drug>1000;

if any(idx0)
    plot( get(gca, 'XLim'), [exinfo.ratemn(idx0) exinfo.ratemn(idx0)], ...
        'Color', c);
end

if any(idx2)
    plot( get(gca, 'XLim'), [exinfo.ratemn_drug(idx2) exinfo.ratemn_drug(idx2)], ...
        'Color', c, 'LineStyle', '--');
end
end