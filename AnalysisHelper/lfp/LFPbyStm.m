function [para] = LFPbyStm(ex)
% compute LFP parameters in the specified ex-file, split by types of
% stimuli
% - ex file must have LFP data, using 'loadCluster.m'
% - ex file must be filtered by 'filterLFP.m'

% stimulus type ==============
[ stimparam, vals] = getStimParam( ex );
% blank_i =  vals > 1;
lenv = length(vals);

% initialization ==================
lfplen = length(ex.Trials(end).LFP_prepro);
lfp_avg = zeros(lenv, lfplen);

fidx = ex.Trials(end).FREQ>= 0 & ex.Trials(end).FREQ<=100;
freq = ex.Trials(end).FREQ(fidx);
lenp = length(ex.Trials(end).POW);
pow_avg = zeros(lenv, lenp);

S = cell(1, lenv);
t = cell(1, lenv);
f = cell(1, lenv);
params = define_params;

C = cell(1, lenv);
phi = cell(1, lenv);
S12 = cell(1, lenv);
S1 = cell(1, lenv);
S2 = cell(1, lenv);
fg = cell(1, lenv);

wnd = 0.064; % was 0.3
wndrange = -wnd:1/1000:wnd;
% ampwindow = wndrange > -0.015 &  wndrange  < 0.015;
ampwindow = wndrange > -wnd &  wndrange  < wnd;
for i = 1:lenv    
    trials = ex.Trials([ex.Trials.(stimparam)] == vals(i));

    % spectrogram
    lfp = vertcat(trials.LFP_prepro);
    lfp = lfp(mean(isnan(lfp), 2)==0, :);
    lfp_cut = lfp(:, 401:end)'; % exclude putative stimulus driven component
    [S{i},t{i},f{i}] = mtspecgramc(lfp_cut, [0.5 0.05], params);
    t{i} = t{i} + ex.Trials(end).LFP_prepro_time(1);
    
    % spike-LFP coherency
    spk = getSpks(trials);
    spk = spk(mean(isnan(lfp),2)==0);
    [C{i},phi{i},S12{i},S1{i},S2{i},fg{i}]= coherencycpt(lfp_cut, cell2struct(spk, 'spk', 1), params);
    
    % LFP averaged trials across the same stimulus
    lfp_avg(i,:) = mean(lfp, 1);

    % delta band
    para.lfp_stm_wave(1).mean(i,:) = nanmean(vertcat(trials.lfp_delta_tc), 1);
    para.lfp_stm_wave(1).pow(i) = nanmean([trials.lfp_delta_pow]);
    
    % theta band
    para.lfp_stm_wave(2).mean(i,:) = nanmean(vertcat(trials.lfp_theta_tc), 1);
    para.lfp_stm_wave(2).pow(i) = nanmean([trials.lfp_theta_pow]);

    % alpha band
    para.lfp_stm_wave(3).mean(i,:) = nanmean(vertcat(trials.lfp_alpha_tc), 1);
    para.lfp_stm_wave(3).pow(i) = nanmean([trials.lfp_alpha_pow]);

    % beta band
    para.lfp_stm_wave(4).mean(i,:) = nanmean(vertcat(trials.lfp_beta_tc), 1);
    para.lfp_stm_wave(4).pow(i) = nanmean([trials.lfp_beta_pow]);

    % gamma band
    para.lfp_stm_wave(5).mean(i,:) = nanmean(vertcat(trials.lfp_gamma_tc), 1);
    para.lfp_stm_wave(5).pow(i) = nanmean([trials.lfp_gamma_pow]);
    
    % averaged power across the same stimulus
    pow_avg(i,:) = nanmean(horzcat(trials.POW),2)';
    
    % spike-triggered LFP
    ex_temp = ex;
    ex_temp.Trials = trials;
     [para.stlfp.avg_stlfp(i,:), para.stlfp.sem_stlfp(i,:), para.stlfp.accspk(i), ...
            para.stlfp.pow(i,:), para.stlfp.freq(i,:), para.stlfp.band(i,:)] = spktriglfp(ex_temp);
%      % CL's correction
%        para.stlfp.avg_stlfp(i,:) = para.stlfp.avg_stlfp(i,:) ...
%            - mean(para.stlfp.avg_stlfp(i, wndrange < -0.06));
       [para.stlfp.peak_stlfp(i), para.stlfp.t_peak_stlfp(i)] =  ...
           getPeakAmplitude(para.stlfp.avg_stlfp(i, ampwindow), wndrange(ampwindow));       
end

% structure
para.stm.param = stimparam;
para.stm.vals = vals;
para.ts = ex.time;
para.ts_cut = ex.time_cut;
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
para.pow_stm.mean = pow_avg(:,fidx);

