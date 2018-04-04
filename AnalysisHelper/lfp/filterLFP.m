function [ex] = filterLFP(ex)
% filter LFP traces into 'band'

% sampling rate
fs = 1000;

% filters
filters = filterSet(fs);

% initial cutoff (0 as no cutoff)
ic = 200;

% assign each trial
for i = 1:length(ex.Trials)
    % delta
    ex.Trials(i).lfp_delta_tc = filtfilt(filters.delta.b, filters.delta.a, ex.Trials(i).LFP_prepro(ic+1:end));
    range = ex.Trials(i).FREQ > 0 & ex.Trials(i).FREQ <4;
    ex.Trials(i).lfp_delta_pow = mean(ex.Trials(i).POW(range));
    
    % theta
    ex.Trials(i).lfp_theta_tc = filtfilt(filters.theta.b, filters.theta.a, ex.Trials(i).LFP_prepro(ic+1:end));
    range = ex.Trials(i).FREQ >=4 & ex.Trials(i).FREQ <=7;
    ex.Trials(i).lfp_theta_pow = mean(ex.Trials(i).POW(range));
    
    % alpha
    ex.Trials(i).lfp_alpha_tc = filtfilt(filters.alpha.b, filters.alpha.a, ex.Trials(i).LFP_prepro(ic+1:end));
    range = ex.Trials(i).FREQ >=8 & ex.Trials(i).FREQ <=13;
    ex.Trials(i).lfp_alpha_pow = mean(ex.Trials(i).POW(range));
    
    % beta
    ex.Trials(i).lfp_beta_tc = filtfilt(filters.beta.b, filters.beta.a, ex.Trials(i).LFP_prepro(ic+1:end));
    range = ex.Trials(i).FREQ >=14 & ex.Trials(i).FREQ <=29;
    ex.Trials(i).lfp_beta_pow = mean(ex.Trials(i).POW(range));
    
    % gamma
    ex.Trials(i).lfp_gamma_tc = filtfilt(filters.gamma.b, filters.gamma.a, ex.Trials(i).LFP_prepro(ic+1:end));
    range = ex.Trials(i).FREQ >=30 & ex.Trials(i).FREQ <=80;
    ex.Trials(i).lfp_gamma_pow = mean(ex.Trials(i).POW(range));
end

% filter for bands
function [filters] = filterSet(fs)
% wave range
delta_range = [0.2, 4];
theta_range = [4,7];
alpha_range = [8, 13];
beta_range = [14, 29];
gamma_range = [30, 80];

% design filter
filters = struct('delta', [], 'theta', [], 'alpha', [], 'beta', [], 'gamma', []);
[filters.delta.b, filters.delta.a] = butter(2, delta_range/(fs/2), 'bandpass');
[filters.theta.b, filters.theta.a] = butter(2, theta_range/(fs/2), 'bandpass');
[filters.alpha.b, filters.alpha.a] = butter(2, alpha_range/(fs/2), 'bandpass');
[filters.beta.b, filters.beta.a] = butter(2, beta_range/(fs/2), 'bandpass');
[filters.gamma.b, filters.gamma.a] = butter(2, gamma_range/(fs/2), 'bandpass');