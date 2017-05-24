function [BRI, ISI_frct] = getBRI_prefStim(exinfo, ex0, ex2)
% computes the burst index as in Anderson et al. 2013
% I only look for the stimuli response that also elicited the highest
% response, i.e. that is closest to the preferred orientation


% declare number of lags (in ms)
lags = 100;

% assign preferred stimuli index and color to plot
idx0 = exinfo.pfi;
idx2 = exinfo.pfi_drug;


%% compute shuffle predictor
% = cross-correlation across all trials with preferred stimulus in 
%  both conditions (5HT/Baseline)

trials = [ex0.raster{idx0}; ex2.raster{idx2}];    
i = 1;
for tr1 = 1:size(trials, 1)
    for tr2 = tr1+1:size(trials, 1)
        xc(:, i) = crosscorr(trials(tr1, :), trials(tr2, :), lags);
        i = i+1;
    end
end

shufflepred_mn = nanmean(xc(lags+1:end, :), 2);
shufflepred_sd = nanstd(xc(lags+1:end, :), 0, 2);

%% get autocorrelation
[acf0, isi0] = getACF(ex0.raster{idx0}, lags, shufflepred_mn, shufflepred_sd);
[acf2, isi2] = getACF(ex2.raster{idx2}, lags, shufflepred_mn, shufflepred_sd);


BRI(1) = sum(nanmean(acf0(2:6, :), 2)) ;
BRI(2) = sum(nanmean(acf2(2:6, :), 2)) ;

ISI_frct(1) = sum(isi0<0.005)/length(isi0); 
ISI_frct(2) = sum(isi2<0.005)/length(isi2); 


%%% plot results
if p_flag
    
    figure('Name',  exinfo.figname);
    col = getCol(exinfo);
    
    subplot(2,1,1); % ISI
    isi0 = isi0(isi0<0.1);
    isi2 = isi2(isi2<0.1);
    h0 = histogram(isi0, 50); hold on;
    h0.FaceColor = 'b';
    h1 = histogram(isi2, h0.BinEdges); hold on;
    h1.FaceColor = col;
    xlabel('ISI [s]'); xlim([0, 0.1]);
    title( sprintf(' base_{isi<5ms}=%1.2f  drug_{isi<5ms}=%1.2f', ISI_frct) );
    
    
    subplot(2,1,2); % ACF
    plot(nanmean(acf0, 2), 'b'); hold on
    plot(nanmean(acf2, 2), col);
    plot(shufflepred_mn);
    plot(shufflepred_mn+shufflepred_sd, ':');
    plot(shufflepred_mn-shufflepred_sd, ':');
    legend('base', exinfo.drugname);
    xlim([1 lags]);
    xlabel('time in ms'); ylabel( 'norm autocorr' );
    title( sprintf(' base_{bri}=%1.2f  drug_{bri}=%1.2f', BRI) );
    
    savefig(gcf, exinfo.fig_bri);
    close(gcf);
    
end

end


function [acf, isi] = getACF(trials, lags, shufflepred_mn, shufflepred_sd)
% returns the autocorrelation for each trial normalized with the shuffle
% predictor, and the interspike intervals that are  across all
% trials

isi = [];

for i = 1:size(trials, 1)
    acf(:, i) = autocorr(trials(i, :), lags);
    acf(:, i) = (acf(:, i) - shufflepred_mn) ./ shufflepred_sd;
    
    % alternatively 
    spks = find(trials(i,:))/1000;
    for j = 1:length(spks)-1
        isi = [isi spks(j+1:end) - spks(j)];
    end
end



end
