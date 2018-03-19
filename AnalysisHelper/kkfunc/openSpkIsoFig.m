function openSpkIsoFig(exinfo)

figdir = 'Z:\Katsuhisa\interaction_project\dataset\Figures\Spike_isolation\';
bad_isolations = [160   163   164   174   202   203   421   562   564   582   819   880];
counts = 0;
for i = 1:length(bad_isolations)
    if exinfo(bad_isolations(i)).isolation_distance(2) < 20 || ...
            exinfo(bad_isolations(i)).Lratio(2) > 0.1
        openfig([figdir exinfo(bad_isolations(i)).figname '_spkiso_base.fig'])
        counts = counts + 1;
    end
    if exinfo(bad_isolations(i)).isolation_distance_drug(2) < 20 || ...
            exinfo(bad_isolations(i)).Lratio_drug(2) > 0.1
        openfig([figdir exinfo(bad_isolations(i)).figname '_spkiso_drug.fig'])
        counts = counts + 1;
    end
end

disp(['the number of experiments with bad isolations: ' num2str(counts)])