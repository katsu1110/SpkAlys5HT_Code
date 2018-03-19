function [exinfo_new] = addspkisolation(exinfo, savefig)

% % load list of index
% load('Z:\Corinna\SharedCode\Katsu\incl_i_one_stimulus_cond.mat')
savedir = 'Z:\Katsuhisa\interaction_project\dataset\Figures\Spike_isolation\';

for i = 1:length(exinfo)
    exinfo(i).snr = [];
    exinfo(i).isi_fraction = [];
    exinfo(i).Lratio = [];
    exinfo(i).isolation_distance = [];

    exinfo(i).snr_drug = [];
    exinfo(i).isi_fraction_drug = [];
    exinfo(i).Lratio_drug = [];
    exinfo(i).isolation_distance_drug = [];
    
    try
        if savefig==0
            [data] = spikeisolation(exinfo(i).fname, 0);
            exinfo(i).snr = data.snr;
            exinfo(i).isi_fraction = data.fraction;
            exinfo(i).Lratio = data.Lratio;
            exinfo(i).isolation_distance = data.isolation_distance;
            
            [data] = spikeisolation(exinfo(i).fname_drug, 0);
            exinfo(i).snr_drug = data.snr;
            exinfo(i).isi_fraction_drug = data.fraction;
            exinfo(i).Lratio_drug = data.Lratio;
            exinfo(i).isolation_distance_drug = data.isolation_distance;
            
            disp(['row ' num2str(i) ' was done!'])
        
        else        
            [data, h] = spikeisolation(exinfo(i).fname, 1);
            exinfo(i).snr = data.snr;
            exinfo(i).isi_fraction = data.fraction;
            exinfo(i).Lratio = data.Lratio;
            exinfo(i).isolation_distance = data.isolation_distance;
            
            savefig(h, [savedir exinfo(i).figname '_spkiso_base.fig'])
            close all;
            
            [data, h] = spikeisolation(exinfo(i).fname_drug, 1);
            exinfo(i).snr_drug = data.snr;
            exinfo(i).isi_fraction_drug = data.fraction;
            exinfo(i).Lratio_drug = data.Lratio;
            exinfo(i).isolation_distance_drug = data.isolation_distance;

            savefig(h, [savedir exinfo(i).figname '_spkiso_drug.fig'])
            close all;
            
            disp(['row ' num2str(i) ' was done (figures are stored) !'])
        end
    catch
        disp(['row ' num2str(i) ' was skipped due to error...'])
    end
    
end
exinfo_new = exinfo;

save('Z:\Katsuhisa\interaction_project\dataset\Data\exinfo_new.mat', 'exinfo_new', '-v7.3'); 