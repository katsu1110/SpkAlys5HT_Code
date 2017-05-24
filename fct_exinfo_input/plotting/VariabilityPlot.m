function  h = VariabilityPlot( exinfo, ex_base, ex_drug )


h = figure('Name', exinfo.figname, 'Position', [509   386   964   596]);
drugcolor = getCol(exinfo);

cax = [-3 3];
kernel = [1 1 1]; kernel = kernel/sum(kernel);

% baseline cluster 0
fnamec0 = strrep(exinfo.fname, exinfo.cluster, 'c0'); % load cluster 0
ex_base_c0 = loadCluster( fnamec0, 'ocul', exinfo.ocul);

[ex_base_c0, spkrate_base_c0] = znormex(ex_base_c0, exinfo); % z norm


% drug cluster 0
fnamec0 = strrep(exinfo.fname_drug, exinfo.cluster, 'c0'); % load cluster 0
ex_drug_c0 = loadCluster( fnamec0, 'ocul', exinfo.ocul);

[ex_drug_c0, spkrate_drug_c0] = znormex(ex_drug_c0, exinfo); % z norm



%%% PLOT
% 1. z-scored data
axes('Position', [0.1 0.6 0.3 0.3]);
scatter([ex_base_c0.Trials.zspkcount], [ex_base.Trials.zspkcount], 20, drugcolor, 'filled', 'MarkerFaceAlpha', 0.7); hold on;
scatter([ex_drug_c0.Trials.zspkcount], [ex_drug.Trials.zspkcount], 20, drugcolor, 'MarkerFaceAlpha', 0.7); hold on;

leg = legend('base', exinfo.drugname, 'Location', 'EastOutside'); leg.Box = 'off';

axis square;
eqax; unity; crossl;

xlabel('z(MU)'); ylabel('z(SU)');
title(sprintf(['Baseline rho=%1.2f, p=%1.2f \n ' exinfo.drugname ' rho=%1.2f, p=%1.2f '], ...
    exinfo.rsc, exinfo.prsc, exinfo.rsc_drug, exinfo.prsc_drug), 'FontSize', 9);


% 2. plot fano facotors
axes('Position', [0.1 0.1 0.3 0.3]);


ff = exinfo.ff.classic;
ff_drug = exinfo.ff_drug.classic;
scatter(ff.spkcnt_mn, ff.spkcnt_var, 20, drugcolor, 'filled', 'MarkerFaceAlpha', 0.7, 'Marker', '^'); hold on;
scatter(ff_drug.spkcnt_mn, ff_drug.spkcnt_var, 20, drugcolor, 'MarkerFaceAlpha', 0.7, 'Marker', '^');

try
    ff2 = exinfo.ff.classic_2ndhalf;
    ff_drug2 = exinfo.ff_drug.classic_2ndhalf;
    scatter(ff2.spkcnt_mn, ff2.spkcnt_var, 20, drugcolor, 'filled', 'MarkerFaceAlpha', 0.7); hold on;
    scatter(ff_drug2.spkcnt_mn, ff_drug2.spkcnt_var, 20, drugcolor, 'MarkerFaceAlpha', 0.7);
catch
    ff2.ff = 0;
    ff_drug2.ff = 0;
end

leg = legend('base all', 'drug all', 'base 2nd half', 'drug 2nd half', ...
    'Location', 'EastOutside'); leg.Box = 'off';

xlabel('spike count mean');
ylabel('spike count variance');

title(sprintf(['all / 2nd half \n Baseline %1.2f / %1.2f, \n ' exinfo.drugname ' %1.2f / %1.2f '], ...
    nanmean(ff.ff), nanmean(ff2.ff), nanmean(ff_drug.ff), nanmean(ff_drug2.ff)), 'FontSize', 9);
set(gca, 'XScale', 'log', 'YScale', 'log')
eqax; unity; crossl; 
axis square





% 3. plot heatmap over time and tuning curves for baseline
colormap('jet');

axes('Position', [0.5 0.8 0.3 0.1]);
imagesc(conv([ex_base_c0.Trials.zspkcount], kernel));     
title('Baseline', 'HorizontalAlignment', 'right');
caxis(cax);    ylabel('MU')
set(gca, 'XTick', [], 'YTick', []);

axes('Position', [0.5 0.65 0.3 0.1]);
imagesc(conv([ex_base.Trials.zspkcount], kernel));  
caxis(cax);    ylabel('SU')
set(gca, 'YTick', []);

c = colorbar('SouthOutside'); c.Position = [0.57    0.5    0.2    0.02];

% tc
axes('Position', [0.85 0.7 0.1 0.1]);
idx = exinfo.ratepar<1000;
plot(exinfo.ratepar(idx), exinfo.ratemn(idx), '+-', 'Color', drugcolor);hold on;
idx = [spkrate_base_c0.(exinfo.param1)]<1000;
plot([spkrate_base_c0(idx).(exinfo.param1)], [spkrate_base_c0(idx).mn], 'o--', 'Color', drugcolor);hold on;
title(sprintf('r_{sig}=%1.2f, p=%1.2f', exinfo.rsig, exinfo.prsig), 'FontSize', 9);

xlabel(exinfo.param1); ylabel('spk/s');



% 4. plot heatmap over time and tuning curves for drug
axes('Position', [0.5 0.25 0.3 0.1]);
imagesc(conv([ex_drug_c0.Trials.zspkcount], kernel));     
title(exinfo.drugname, 'HorizontalAlignment', 'right');
caxis(cax);    ylabel('MU')
set(gca, 'XTick', [], 'YTick', []);

axes('Position', [0.5 0.1 0.3 0.1]);
imagesc(conv([ex_drug.Trials.zspkcount], kernel));

caxis(cax);    set(gca, 'YTick', []);
ylabel('SU'); xlabel('trials');


%tc
axes('Position', [0.85 0.2 0.1 0.1]);
idx = exinfo.ratepar_drug<1000;
plot(exinfo.ratepar_drug(idx), exinfo.ratemn_drug(idx), '+-', 'Color', drugcolor);hold on;
idx = [spkrate_drug_c0.(exinfo.param1)]<1000;
plot([spkrate_drug_c0(idx).(exinfo.param1)], [spkrate_drug_c0(idx).mn], 'o--', 'Color', drugcolor);hold on;
title(sprintf('r_{sig}=%1.2f, p=%1.2f', exinfo.rsig, exinfo.prsig), 'FontSize', 9);

xlabel(exinfo.param1); ylabel('spk/s');
leg = legend('SU', 'MU', 'Orientation', 'horizontal');
leg.Box = 'off';
leg.FontSize = 7;
leg.Position = [0.85,0.08,0.1,0.01];
box off;


% 5. save figure
savefig(h, exinfo.fig_noisecorr);
delete(h);

end


