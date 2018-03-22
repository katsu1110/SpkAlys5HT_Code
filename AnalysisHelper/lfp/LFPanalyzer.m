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

% spike-triggered LFP
[avg_stlfp, sem_stlfp, accspk, pow, freq] = spktriglfp(ex0, 'time', 0.064);
para.cond(1).stlfp = struct('mean', avg_stlfp, 'sem', sem_stlfp, ...
    'acc', accspk, 'pow', pow, 'freq', freq);
[avg_stlfp, sem_stlfp, accspk, pow, freq] = spktriglfp(ex2, 'time', 0.064);
para.cond(2).stlfp = struct('mean', avg_stlfp, 'sem', sem_stlfp, ...
    'acc', accspk, 'pow', pow, 'freq', freq);

% LFP separated by stimulus type
function [para] = LFPbyStm(ex)
% stimulus type ==============
[ stimparam, vals] = getStimParam( ex );
lenv = length(vals);

% aligned time ==================
ts = ex.time;

% initialization ==================
lfplen = length(ex.Trials(end).LFP_prepro);
lfp_avg = zeros(lenv, lfplen);
delta_avg = zeros(lenv, lfplen);
theta_avg = zeros(lenv, lfplen);
alpha_avg = zeros(lenv, lfplen);
beta_avg = zeros(lenv, lfplen);
gamma_avg = zeros(lenv, lfplen);
delta_pow = zeros(lenv, 1);
theta_pow = zeros(lenv, 1);
alpha_pow = zeros(lenv, 1);
beta_pow = zeros(lenv, 1);
gamma_pow = zeros(lenv, 1);

fidx = ex.Trials(end).FREQ>= 0 & ex.Trials(end).FREQ<=100;
freq = ex.Trials(end).FREQ(fidx);
lenp = length(ex.Trials(end).POW);
pow_avg = zeros(lenv, lenp);

S = cell(1, lenv);
t = cell(1, lenv);
f = cell(1, lenv);
params.err = 0;
params.Fs = 1000;
params.fpass = [0 100];
params.tapers = [2,3];
params.pad = 0;
params.trialave = 1;

C = cell(1, lenv);
phi = cell(1, lenv);
S12 = cell(1, lenv);
S1 = cell(1, lenv);
S2 = cell(1, lenv);
fg = cell(1, lenv);

for i = 1:lenv    
    trials = ex.Trials([ex.Trials.(stimparam)] == vals(i));

    % spectrogram
    lfp = vertcat(trials.LFP_prepro);
    lfp = lfp(mean(isnan(lfp), 2)==0, :);
    [S{i},t{i},f{i}] = mtspecgramc(lfp', [0.5 0.05], params);
    t{i} = t{i} + ex.Trials(end).LFP_prepro_time(1);
    
    % spike-LFP coherency
    spk = getSpks(trials);
    spk = spk(mean(isnan(lfp),2)==0);
    [C{i},phi{i},S12{i},S1{i},S2{i},fg{i}]= coherencycpt(lfp',spk,params);
    
    % LFP averaged trials across the same stimulus
    lfp_avg(i,:) = mean(lfp, 1);

    % delta band
    delta_avg(i,:) = nanmean(vertcat(trials.lfp_delta_tc), 1);
    delta_pow(i) = nanmean([trials.lfp_delta_pow]);
    
    % theta band
    theta_avg(i,:) = nanmean(vertcat(trials.lfp_theta_tc), 1);
    theta_pow(i) = nanmean([trials.lfp_theta_pow]);

    % alpha band
    alpha_avg(i,:) = nanmean(vertcat(trials.lfp_alpha_tc), 1);
    alpha_pow(i) = nanmean([trials.lfp_alpha_pow]);

    % beta band
    beta_avg(i,:) = nanmean(vertcat(trials.lfp_beta_tc), 1);
    beta_pow(i) = nanmean([trials.lfp_beta_pow]);

    % gamma band
    gamma_avg(i,:) = nanmean(vertcat(trials.lfp_gamma_tc), 1);
    gamma_pow(i) = nanmean([trials.lfp_gamma_pow]);
    
    % averaged power across the same stimulus
    pow_avg(i,:) = nanmean(horzcat(trials.POW),2)';
    
    % spike-triggered LFP
    ex_temp = ex;
    ex_temp.Trials = trials;
    [para.stlfp.avg_stlfp(i,:), para.stlfp.sem_stlfp(i,:), para.stlfp.accspk(i), ...
        para.stlfp.pow(i,:), para.stlfp.freq(i,:), para.stlfp.band(i,:)] = spktriglfp(ex_temp, 'time', 0.064);
end

% structure
para.stm.param = stimparam;
para.stm.vals = vals;
para.ts = ts;
para.f = freq';
para.params = params;
para.spectrogram.S = S;
para.spectrogram.t = t;
para.spectrogram.f = f;
para.coherence.C = C;
para.coherence.phi = phi;
para.coherence.S12 = S12;
para.coherence.S1 = S1;
para.coherence.S2 = S2;
para.coherence.f = fg;
para.lfp_stm.mean = lfp_avg;
para.lfp_stm_wave(1).mean = delta_avg;
para.lfp_stm_wave(2).mean = theta_avg;
para.lfp_stm_wave(3).mean = alpha_avg;
para.lfp_stm_wave(4).mean = beta_avg;
para.lfp_stm_wave(5).mean = gamma_avg;
para.lfp_stm_wave(1).power = delta_pow;
para.lfp_stm_wave(2).power = theta_pow;
para.lfp_stm_wave(3).power = alpha_pow;
para.lfp_stm_wave(4).power = beta_pow;
para.lfp_stm_wave(5).power = gamma_pow;
para.pow_stm.mean = pow_avg(:,fidx);

% get spikes given trials
function spk_out = getSpks(trials)

for i = 1:length(trials)
    % spikes within the stimulus presentation time
    t_strt = trials(i).Start - trials(i).TrialStart;
    t_end = t_strt(end)+mean(diff(t_strt));

    spk = trials(i).Spikes( trials(i).Spikes >= t_strt(1) & ...
        trials(i).Spikes <= t_end) - t_strt(1); 
    spk = round(spk*1000)/1000;

    spk_out(i).spk = spk;
end

% visualize
function visualizer(para, name0, name2, prefix, save_flag)
name = {name0, name2};
savedir = 'Z:\Katsuhisa\serotonin_project\LFP_project\Figures_kk\';

% spike-triggered LFP =========================
h = figure;
wnd = 0.064; 
for k = 1:2
    switch k
        case 1
            col = zeros(1,3);
        case 2
            col = [1,0,0];
    end
    
    % stlfp
    subplot(1,2,1)
    fill_between(-wnd:0.001:wnd, para.cond(k).stlfp.mean - para.cond(k).stlfp.sem, ...
        para.cond(k).stlfp.mean + para.cond(k).stlfp.sem, col);
    hold on;
    plot(-wnd:0.001:wnd, para.cond(k).stlfp.mean, '-','color',col);
    hold on;
    
    % power and frequency
    subplot(1,2,2)
    plot(para.cond(k).stlfp.freq, para.cond(k).stlfp.pow, '-','color',col);
    hold on;
end

% cosmetics
subplot(1,2,1)
yy = get(gca, 'YLim');
hold on;
plot([0 0], yy, '-k')
xlabel('time rel:spike [s]');
ylabel('avg LFP +/- SEM (\muV)');
xlim([-wnd wnd]);
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

subplot(1,2,2)
xlabel('frequency [Hz]');
ylabel('power');
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
xlim([0 100])

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
            plot(para.cond(k).lfpstm.ts, para.cond(k).lfpstm.lfp_stm_wave(b).mean(s,:), ...
                marker, 'color', col(s,:))
            xlim([-0.2 max(para.cond(k).lfpstm.ts(~isnan(para.cond(k).lfpstm.ts)))])
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
        para.cond(1).lfpstm.lfp_stm_wave(b).power(para.cond(1).lfpstm.stm.vals < 1000), ...
        '-or')
    hold on;
    plot(para.cond(2).lfpstm.stm.vals(para.cond(2).lfpstm.stm.vals < 1000), ...
        para.cond(2).lfpstm.lfp_stm_wave(b).power(para.cond(2).lfpstm.stm.vals < 1000), ...
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
        polarhistogram(para.cond(1).lfpstm.coherence.phi{i}, 'FaceColor','red','FaceAlpha',.3);
        title([para.cond(1).lfpstm.stm.param, ' = ' ...
               num2str(para.cond(1).lfpstm.stm.vals(i))])
    catch
        disp(['stimulus index ' num2str(i) ' error'])
    end
    try
        subplot(2, lenv, i+lenv)
        polarhistogram(para.cond(2).lfpstm.coherence.phi{i}, 'FaceColor','red','FaceAlpha',.3);
    catch
        disp(['stimulus index ' num2str(i) ' error'])
    end
end

figname = [name{1}, 'VS', name{2}, '_' 'LfpPhaseOfSpikes_' prefix];
set(h, 'Name', figname,'NumberTitle','off')
if save_flag==1
    savefig(h, strcat(savedir, figname, '.fig'))
end
