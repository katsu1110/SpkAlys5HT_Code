function spectogramPlot(exinfo)
%spectogramPlot returns a plot with a time-frequency spectogram for each
%condition and their difference

% addpath(genpath('C:\Users\Corinna\Documents\chronux_2_11'), '-end');

lfpname = strrep(exinfo.fname, exinfo.cluster, 'lfp');
lfpname = strrep(lfpname, '.mat', [num2str(exinfo.ocul) '.mat']);
ex0 = loadCluster(lfpname);


lfpname2 = strrep(exinfo.fname_drug, exinfo.cluster, 'lfp');
lfpname2 = strrep(lfpname2, '.mat', [num2str(exinfo.ocul) '.mat']);
ex2 = loadCluster(lfpname2);

if ~isempty(ex0.Trials) & ~isempty(ex2.Trials)
    
    
    ind = [ex0.Trials.me] == exinfo.ocul & ...
        [ex0.Trials.(exinfo.param1)] == exinfo.ratepar(exinfo.pfi);
    ppLFP1 = vertcat(ex0.Trials(ind).LFP_interp);
    
    
    ind2 = [ex2.Trials.me] == exinfo.ocul & ...
        [ex2.Trials.(exinfo.param1)] == exinfo.ratepar_drug(exinfo.pfi_drug);
    ppLFP2 = vertcat(ex2.Trials(ind2).LFP_interp);
    
    
    [S1,T1,F1] = spectogramPlot_Helper( ppLFP1, 1000 );
    [S2,T2,F2] = spectogramPlot_Helper( ppLFP2, 1000 );
    
    
    
    %%
    %4.) Plot Power Specturm and LFP
    h = figure('Name', exinfo.figname, 'Position', [351 531 1342 420] );
    
    time = 0:0.001:0.5;
    c = getCol(exinfo);
    
    %%% Baseline
    subplot(2, 4, 1);    % mean and std lfp
    plot(time, nanmean(ppLFP1), 'b', time, ...
        nanmean(ppLFP1)+ nanstd(ppLFP1),  'b:',....
        time, nanmean(ppLFP1)- nanstd(ppLFP1), 'b:', ...
        time, nanmean(ppLFP2), c); ho
    plot(time, nanmean(ppLFP2)+ nanstd(ppLFP2), ':', 'Color', c);
    plot(time, nanmean(ppLFP2)- nanstd(ppLFP2), ':', 'Color', c);
    xlim([min(time) max(time)]);
    xlabel('time [s]'); ylabel('\muVolt');
    title(['baseline(b) vs ' exinfo.drugname '(' c ')']);
    
    
    subplot(2, 4, 5);    % power spect
    plot(F1, mean(mean(S1,1), 3), 'b'); hold on;
    plot(F2, mean(mean(S2,1), 3), c);
    xlabel('frequency'); ylabel('spectrum');
    
    subplot(2, 4, [2 6]);   % time-spectogram baseline
    imagesc(T1, F1, mean(S1, 3)'); colorbar;
    ylim([0 100]); xlim([min(time) max(time)]);
    title(sprintf('Baseline LFP spectogram (n= %1d)', size(ppLFP1,2)) );
    xlabel('time'); ylabel('frequency');
    
    subplot(2, 4, [3 7]);   % time-spectogram drug
    imagesc(T2, F2, mean(S2, 3)'); colorbar;
    ylim([0 100]); xlim([min(time) max(time)]);
    title(sprintf([exinfo.drugname ' LFP spectogram (n= %1d)'], size(ppLFP2,2)));
    xlabel('time');
    
    
    subplot(2, 4, [4 8]);   % time-spectogram drug
    imagesc(T2, F2, mean(S1, 3)' - mean(S2, 3)'); colorbar;
    ylim([0 100]); xlim([min(time) max(time)]);
    title('LFP spectogram difference');
    xlabel('time');
    
    savefig(h, exinfo.fig_lfpSpec);
    close(h);
else
    
end
end




function [S,T,F] = spectogramPlot_Helper( rawEphys, samplingRate )


addpath( genpath('C:\Users\Corinna\Documents\chronux_2_11'),'-end')

%1.) Define Variables
time = 0:0.001:0.5;
targetSamplingRate = samplingRate;
% samplingRateReduction = round(samplingRate/targetSamplingRate);
WindowSize = 0.05;
WindowStep = WindowSize/3;
params.Fs = targetSamplingRate;
params.fpass = [0 100];
params.trialave = 1;
params.tapers = [7 6];

%2.) Get LFP
%Decimate Trace
% rawLFP = decimate(rawEphys,samplingRateReduction);
%Remove linear trends
% rawLFP = detrend(rawLFP);
rawLFP = rawEphys;

%3.) Get Spectrogram
[S,T,F] = mtspecgramc(rawLFP',[WindowSize,WindowStep],params);
T = T + min(time);

rmpath(genpath('C:\Users\Corinna\Documents\chronux_2_11'));

end


