function [para] = LFPanalyzer(exinfo, varargin)
%% comprehensive LFP analysis
% INPUT: exinfo (but assumes a pair of sessions)
%
% OUTPUT: lfpdata ... structure containing analized data
% 
% written by Katsuhisa (19.03.18)
% +++++++++++++++++++++++++++++++++++++++

% deal with input arguments
plot_flag = 0;
save_flag = 0;
j = 1;
while j <= nargin - 1
    switch varargin{j}
        case 'plot'
            plot_flag = 1;
            j = j + 1;
        case 'save'
            save_flag = 1;
            j = j + 1;
    end
end
    
% stimulus & drugname
para.stimulus = exinfo.param1;
if exinfo.isRC
    para.stimulus = 'rc';
end
para.drugname = exinfo.drugname;

% load lfp data and preprocessing 
% (POW, FREQ, LFP_prepro, LFP_prepro_time)
ex0 = loadCluster(exinfo.fname, 'loadlfp',1);
ex2 = loadCluster(exinfo.fname_drug, 'loadlfp',1);

% median split of the pupil size
[ex0_sps, ex0_lps, ex0] = pupilSplit(ex0);
[ex2_sps, ex2_lps, ex2] = pupilSplit(ex2);

% LFP analysis
[para.drug] = LFPanalyzerWrapper(ex0, ex2);
[para.ps_base] = LFPanalyzerWrapper(ex0_sps, ex0_lps);
[para.ps_drug] = LFPanalyzerWrapper(ex2_sps, ex2_lps); 

% visualization
if plot_flag==1
    close all;
    visualizer(para.drug, 'base', para.drugname, exinfo.figname, save_flag)
    visualizer(para.ps_base, 'S-ps_base', 'L-ps_base', exinfo.figname, save_flag)
    visualizer(para.ps_drug, ['S-ps_' para.drugname], ['L-ps_' para.drugname], exinfo.figname, save_flag)
end


%% subfunctions
% wrapper
function [para] = LFPanalyzerWrapper(ex0, ex2)
% filter LFP traces into bands
[ex0] = filterLFP(ex0);
[ex2] = filterLFP(ex2);

% LFP properties across stimulus type
[para.cond(1).lfpstm] = LFPbyStm(ex0);
[para.cond(2).lfpstm] = LFPbyStm(ex2);

% % spike-triggered LFP --- already done in each stimulus type
% [avg_stlfp, sem_stlfp, accspk, pow, freq] = spktriglfp(ex0);
% para.cond(1).stlfp = struct('mean', avg_stlfp, 'sem', sem_stlfp, ...
%     'acc', accspk, 'pow', pow, 'freq', freq);
% [avg_stlfp, sem_stlfp, accspk, pow, freq] = spktriglfp(ex2);
% para.cond(2).stlfp = struct('mean', avg_stlfp, 'sem', sem_stlfp, ...
%     'acc', accspk, 'pow', pow, 'freq', freq);


% visualize
function visualizer(para, name0, name2, prefix, save_flag)
name = {name0, name2};
if mean(ismember('gpfs0', cd))==1
    savedir = '/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/Figures_kk/';
else
    savedir = 'Z:\Katsuhisa\serotonin_project\LFP_project\Figures_kk\';
end

% spike-triggered LFP =========================
h = figure;
wnd = 0.1; % was 0.3
for k = 1:2
    switch k
        case 1
            col = zeros(1,3);
        case 2
            col = [1,0,0];
    end
    
    % stlfp
    subplot(1,3,1)
    fill_between(-wnd:0.001:wnd, para.cond(k).lfpstm.stlfp.avg_stlfp(end,:) - ...
        para.cond(k).lfpstm.stlfp.sem_stlfp(end,:), ...
        para.cond(k).lfpstm.stlfp.avg_stlfp(end,:) + para.cond(k).lfpstm.stlfp.sem_stlfp(end,:), col);
    hold on;
    plot(-wnd:0.001:wnd, para.cond(k).lfpstm.stlfp.avg_stlfp(end,:), '-','color',col);
    hold on;
    
    % power and frequency (< 30Hz)
    freq = para.cond(k).lfpstm.stlfp.freq(end,:);
    subplot(1,3,2)
    plot(freq(freq<30), para.cond(k).lfpstm.stlfp.pow(end,freq<30), '-','color',col);
    hold on;
    
    % power and frequency (>= 30Hz)
    subplot(1,3,3)
    plot(freq(freq>=30), para.cond(k).lfpstm.stlfp.pow(end,freq>=30), '-','color',col);
    hold on;
end

% cosmetics
subplot(1,3,1)
yy = get(gca, 'YLim');
hold on;
plot([0 0], yy, '-k')
xlabel('time rel:spike [s]');
ylabel('avg LFP +/- SEM (\muV)');
xlim([-wnd wnd]);
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

subplot(1,3,2)
xlabel('frequency [Hz]');
ylabel('power');
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
xlim([0 30])

subplot(1,3,3)
xlabel('frequency [Hz]');
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
xlim([30 100])

% stlfp0 = para.cond(1).stlfp.mean([-wnd:0.001:wnd]==0);
% stlfp2 = para.cond(2).stlfp.mean([-wnd:0.001:wnd]==0);

figname = [name{1}, 'VS', name{2}, '_' 'spLFP_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end

% stimulus type =====================================
h = figure;    
for k = 1:2       
    lenv = length(para.cond(k).lfpstm.stm.vals);
    col = lines(lenv);
    for s = 1:lenv
        % LFP traces
        subplot(2,2,1+2*(k-1))
%         fill_between(para.cond(k).lfpstm.ts, ...
%             para.cond(k).lfpstm.lfp_stm.mean(s,:) - para.cond(k).lfpstm.lfp_stm.sem(s,:),...
%             para.cond(k).lfpstm.lfp_stm.mean(s,:) + para.cond(k).lfpstm.lfp_stm.sem(s,:),...
%             col(s,:))
%         hold on;
        plot(para.cond(k).lfpstm.ts, para.cond(k).lfpstm.lfp_stm.mean(s,:), '-', 'color', col(s,:))
        hold on;

        % power
        subplot(2,2,2+2*(k-1))
        plot(para.cond(k).lfpstm.f, para.cond(k).lfpstm.pow_stm.mean(s,:), '-', 'color', col(s,:))
        hold on;
    end

    % cosmetics
    subplot(2,2,1+2*(k-1))
    yy = get(gca, 'YLim');
    plot([0 0], yy, '-k')
    if k==2
        xlabel('time (s)')
    end
    message = sprintf([name{k} '\n LFP']);
    xlim([-0.2 max(para.cond(k).lfpstm.ts(~isnan(para.cond(k).lfpstm.ts)))])
    ylabel(message)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    subplot(2,2,2+2*(k-1))
    yy = get(gca, 'YLim');
    plot([0 0], yy, '-k')
    if k==2
        xlabel('frequency (Hz)')
    end
    ylabel('power')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

figname = [name{1}, 'VS', name{2}, '_' 'LFPandPOW_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end

% stimulus tuning by LFP bands ============
h = figure;    
for k = 1:2       
    switch k
        case 1
            marker = '-';
        case 2
            marker = '--';
    end
    lenv = length(para.cond(k).lfpstm.stm.vals(para.cond(k).lfpstm.stm.vals < 1000));
    col = lines(lenv);
    for s = 1:lenv
        if para.cond(k).lfpstm.stm.vals(s) >= 1000
            continue
        end
        for b = 1:5
            % LFP traces
            subplot(2,5,b)
            plot(para.cond(k).lfpstm.ts_cut, para.cond(k).lfpstm.lfp_stm_wave(b).mean(s,:), ...
                marker, 'color', col(s,:))
            xlim([-0.2 max(para.cond(k).lfpstm.ts_cut(~isnan(para.cond(k).lfpstm.ts_cut)))])
            hold on;
        end
    end
end
for b = 1:5
    switch b
        case 1
            band = 'delta';
        case 2
            band = 'theta';
        case 3
            band = 'alpha';
        case 4
            band = 'beta';
        case 5
            band = 'gamma';
    end
    subplot(2,5,b)
    yy = get(gca, 'YLim');
    plot([0 0], yy, '-k')
    title(band)
    if b==1
        xlabel('time (s)')
        ylabel('LFP')
    end
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    % tuning curve
    subplot(2,5,b+5)
    plot(para.cond(1).lfpstm.stm.vals(para.cond(1).lfpstm.stm.vals < 1000), ...
        para.cond(1).lfpstm.lfp_stm_wave(b).pow(para.cond(1).lfpstm.stm.vals < 1000), ...
        '-or')
    hold on;
    plot(para.cond(2).lfpstm.stm.vals(para.cond(2).lfpstm.stm.vals < 1000), ...
        para.cond(2).lfpstm.lfp_stm_wave(b).pow(para.cond(2).lfpstm.stm.vals < 1000), ...
        '--^r')
    if b==1
        ylabel('power')
        xlabel(para.cond(1).lfpstm.stm.param)
    end
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

figname = [name{1}, 'VS', name{2}, '_' 'LFPtuning_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end


% stimulus tuning by stLFP bands ============
h = figure;    

for b = 1:5
    switch b
        case 1
            band = 'delta';
        case 2
            band = 'theta';
        case 3
            band = 'alpha';
        case 4
            band = 'beta';
        case 5
            band = 'gamma';
    end
    subplot(1,5,b)
    yy = get(gca, 'YLim');
    plot([0 0], yy, '-k')    
    
    % tuning curve
    subplot(1,5,b)
    plot(para.cond(1).lfpstm.stm.vals(para.cond(1).lfpstm.stm.vals < 1000), ...
        para.cond(1).lfpstm.stlfp.band(para.cond(1).lfpstm.stm.vals < 1000, b), ...
        '-or')
    hold on;
    plot(para.cond(2).lfpstm.stm.vals(para.cond(2).lfpstm.stm.vals < 1000), ...
        para.cond(2).lfpstm.stlfp.band(para.cond(2).lfpstm.stm.vals < 1000, b), ...
        '--^r')
    if b==1
        ylabel('power')
        xlabel(para.cond(1).lfpstm.stm.param)
        title(band)
    end
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

figname = [name{1}, 'VS', name{2}, '_' 'stLFPtuning_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end


% Spectrogram ===================== 
h = figure;
crange = [0, 1];
for k = 1:2        
    for s = 1:lenv        
        subplot(3, lenv, s + lenv*(k-1))
        imagesc(para.cond(k).lfpstm.spectrogram.t{s}, ...
            para.cond(k).lfpstm.spectrogram.f{s},...
            mag2db(para.cond(k).lfpstm.spectrogram.S{s})')
        axis xy;
        colormap('jet')
        c = caxis;
        colorbar('off')            
        if crange(1) > c(1)
            crange(1) = c(1);
        end
        if crange(2) < c(2)
            crange(2) = c(2);
        end           
    end

    for s = 1:lenv
        subplot(3, lenv, s + lenv*(k-1))
        caxis(crange)
        if k==1
            title([para.cond(k).lfpstm.stm.param, ' = ' ...
                num2str(para.cond(k).lfpstm.stm.vals(s))])
        else
            title('')
        end        
        if s==1
            message = sprintf([name{k} ' \n frequency (Hz)']);
            ylabel(message);
        else
            ylabel('')
        end
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    end
end

% difference
crange = [-1 1];
for s = 1:lenv
    subplot(3, lenv, s + lenv*2)
    imagesc(para.cond(1).lfpstm.spectrogram.t{s}, ...
            para.cond(2).lfpstm.spectrogram.f{s},...
            (mag2db(para.cond(1).lfpstm.spectrogram.S{s})...
            - mag2db(para.cond(2).lfpstm.spectrogram.S{s}))')
    axis xy;
    colormap('jet')
    colorbar('off')
    
    c = caxis;
    if crange(1) > c(1)
        crange(1) = c(1);
    end
    if crange(2) < c(2)
        crange(2) = c(2);
    end

    title('')
    if s==1
        xlabel('time [s]');
        message = sprintf([name{1} ' - ' name{2} ' \n frequency (Hz)']);
        ylabel(message);
    else
        xlabel('')
        ylabel('')
    end
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end
for s = 1:lenv
    subplot(3, lenv, s + lenv*2)
    caxis(crange)
end
figname = [name{1}, 'VS', name{2}, '_' 'Spectrogram_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end

 % Spike-LFP coherency ===================== 
h = figure;
yrange1 = [0 0];
yrange2 = [0 0];
for s = 1:lenv
    subplot(2, lenv, s)
    plot(para.cond(1).lfpstm.coherence.f{s}, para.cond(1).lfpstm.coherence.C{s}, '-', 'color', col(s,:))
    hold on;
    plot(para.cond(2).lfpstm.coherence.f{s}, para.cond(2).lfpstm.coherence.C{s}, '--', 'color', col(s,:))
    title([para.cond(1).lfpstm.stm.param, ' = ' ...
            num2str(para.cond(1).lfpstm.stm.vals(s))])
    yy = get(gca, 'YLim');
    if yrange1(1) > yy(1)
        yrange1(1) = yy(1);
    end
    if yrange1(2) < yy(2)
        yrange1(2) = yy(2);
    end
    if s==1
        ylabel('coherence');
    end
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    subplot(2, lenv, s + lenv)
    plot(para.cond(1).lfpstm.coherence.f{s}, ...
        zeros(1, length(para.cond(1).lfpstm.coherence.f{s})), '-k')
    hold on;
    plot(para.cond(1).lfpstm.coherence.f{s}, para.cond(1).lfpstm.coherence.C{s} - ...
        para.cond(2).lfpstm.coherence.C{s}, ':', 'color', col(s,:))
    yy = get(gca, 'YLim');
    if yrange2(1) > yy(1)
        yrange2(1) = yy(1);
    end
    if yrange2(2) < yy(2)
        yrange2(2) = yy(2);
    end
    if s==1
        xlabel('frequency (Hz)');
        ylabel('\Delta coherence')
    end
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end
for s = 1:lenv
    subplot(2, lenv, s)
    ylim(yrange1)
    
    subplot(2, lenv, s + lenv)
    ylim(yrange2)
end

figname = [name{1}, 'VS', name{2}, '_' 'SpkLfpCoherence_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end

 % Spike-LFP phase ===================== 
h = figure;
for i = 1:lenv
    try
        subplot(2, lenv, i)
        polarhistogram(para.cond(1).lfpstm.coherence.phi{i}, ...
            length(para.cond(1).lfpstm.coherence.phi{i}), ...
            'FaceColor','red','FaceAlpha',.3);
        title([para.cond(1).lfpstm.stm.param, ' = ' ...
               num2str(para.cond(1).lfpstm.stm.vals(i))])
    catch
        disp(['stimulus index ' num2str(i) ' error'])
    end
    try
        subplot(2, lenv, i+lenv)
        polarhistogram(para.cond(2).lfpstm.coherence.phi{i}, ...
            length(para.cond(2).lfpstm.coherence.phi{i}), ...
            'FaceColor','red','FaceAlpha',.3);
    catch
        disp(['stimulus index ' num2str(i) ' error'])
    end
end

figname = [name{1}, 'VS', name{2}, '_' 'LfpPhaseOfSpikes_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end
