function exinfo = getISI_All(exinfo, ex0, ex2, p_flag)
% computes the burst index as in Anderson et al. 2013
% I compute the isi across all stimuli condition.. 



if exinfo.isadapt
    return
else
    time = 0:0.001:0.45;
end


trials0 = ex0.Trials( [ex0.Trials.(exinfo.param1)] == exinfo.ratepar(exinfo.pfi) );
trials2 = ex2.Trials( [ex2.Trials.(exinfo.param1)] == exinfo.ratepar(exinfo.pfi) );

%%% merge inter spike interval from all trials
% baseline
[isi0, acf0, ccf_alltrials0] = getISI(trials0, time);
ccf0 = nanmean(ccf_alltrials0(:,451:end));
sd = nanstd(ccf_alltrials0(:,451:end), 0, 1);         %# after Anderson et al. (2013)
acf_norm0 = (acf0-ccf0) ./ sd;

% drug
[isi2, acf2, ccf_alltrials2]= getISI(trials2, time);
ccf2 = nanmean(ccf_alltrials2(:,451:end));
sd = nanstd(ccf_alltrials2(:,451:end), 0, 1);         
acf_norm2 = (acf2-ccf2) ./ sd;

 
%% assign results
exinfo.isi_frct(1) = sum(isi0<0.005)/length(isi0); 
exinfo.isi_frct(2) = sum(isi2<0.005)/length(isi2); 

exinfo.bridx(1) = mean(acf_norm0(1:5));  % average < 5ms
exinfo.bridx(2) = mean(acf_norm2(1:5)); 



%% plot results
if p_flag
    
    h = figure('Name',  exinfo.figname);
    col = getCol(exinfo);
    
    subplot(3,1,1);
    isi0 = isi0(isi0<0.05);
    isi2 = isi2(isi2<0.05);
    h0 = histogram(isi0*1000, 0:50); hold on;             h0.FaceColor = 'b';
    h1 = histogram(isi2*1000, h0.BinEdges); hold on;     h1.FaceColor = col;
    hold on;    plot([5 5], get(gca, 'YLim'), 'k');

    xlabel('ISI [ms]'); xlim([0, 50]); legend('baseline', exinfo.drugname);
    title( sprintf(['base_{isi<5ms}: %1.1f %%   ' exinfo.drugname '_{isi<5ms}: %1.1f %%'], ...
        exinfo.isi_frct(1)*100, exinfo.isi_frct(2)*100) );

      
    
    
    subplot(3,1,2);
    plot(time*1000, acf0, 'color', 'b'); hold on
    plot(time*1000, acf2, 'color', col);
        
    ylabel('ACF'); xlabel('time [ms]'); xlim([1, 450]); legend('baseline', exinfo.drugname);
    set(gca, 'XScale', 'log', 'XTick', [5 10 20 50 100 200 400])
    
    title( sprintf(['BRI (mean acf 1-4ms): base: %1.1f ' exinfo.drugname ': %1.1f'], ...
        exinfo.bridx(1), exinfo.bridx(2)) );
        
    
    subplot(3,1,3);
    plot(time*1000, acf_norm0, 'color', 'b'); hold on
    plot(time*1000, acf_norm2, 'color', col);
    ylabel('ACF norm [sd(ccf)]'); xlabel('time [ms]'); xlim([1, 10]); 
    legend('baseline', exinfo.drugname);
    set(gca, 'XScale', 'log', 'XTick', [5 10 20 50 100 200 400])
       
    savefig(h, exinfo.fig_bri);
    close(h);
    
end

end


function [isi, acf, ccf] = getISI(trials, time)
% returns the interstimulus intervals, the autocorrelation, and the shuffle
% predictor, that is the cross-correlation across trials
% see Anderson et al. (2013) or Compte et al. (2003)

acf_2 = []; ccf = nan(1, 901) ; isi = [];
raster = spikes2raster(trials, time);


for i = 1:length(trials)

    % select spikes within stimulus presentation range
    t_strt = trials(i).Start - trials(i).TrialStart; % stimuli frame onsets
    ind = trials(i).Spikes > t_strt(1) & trials(i).Spikes < t_strt(end)+0.1;
    spks = trials(i).Spikes(ind);
    

    isi = [isi; diff(spks)];
    
    for j = 1:length(spks)-1
        acf_2 = [acf_2; spks(j+1:end) - spks(j)];
    end
    
    for k = i+1:length(trials)
        ccf = [ccf; xcorr(raster(i, :), raster(k, :))];
    end
    
end

h = histogram(acf_2, 0:0.001:0.451);
acf = h.Values ./ length(trials);  close all;


end

function raster = spikes2raster(trials, time)

raster = zeros(length(trials), length(time));
for t = 1:length(trials)
    t_strt = trials(t).Start - trials(t).TrialStart;

    spk =  trials(t).Spikes(...
        trials(t).Spikes>=t_strt(1) & ...
        trials(t).Spikes<=t_strt(end))-t_strt(1)  ;
    
    % adding it to the raster matrix
    idx = round(( spk + abs(time(1)) ) *1000);
    idx(idx==0) = 1; % avoid bad indexing
    raster(t, idx) = 1;
end


end




