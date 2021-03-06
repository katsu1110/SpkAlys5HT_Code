function [ex] = filterLFP(ex)
% filter LFP traces into 'band'

% sampling rate
fs = 1000;

% filters
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
range = {[0.4,3.9], [4, 7], [8, 13], [14, 29], [30, 80]};
lenb = length(bands);
filtdim = 2;
for i = 1:lenb
    [filters.(bands{i}).b, filters.(bands{i}).a] = butter(filtdim, range{i}/(fs/2), 'bandpass');
end

% assign each trial
for i = 1:length(ex.Trials)    
    for b = 1:lenb
        % all
        ex.Trials(i).(['lfp_' bands{b} '_tc']) = ...
                filtfilt(filters.(bands{b}).b, filters.(bands{b}).a, ex.Trials(i).LFP_z);
        freq = ex.Trials(i).FREQ;
        ex.Trials(i).(['lfp_' bands{b} '_pow']) = ...
            mean(ex.Trials(i).POW(freq >= range{b}(1) & freq <= range{b}(2)));
        % baseline, stimulus evoked, sustained
        for u = 1:3
            ex.Trials(i).period(u).(['lfp_' bands{b} '_tc']) = ...
                filtfilt(filters.(bands{b}).b, filters.(bands{b}).a, ex.Trials(i).period(u).LFP_z);
            freq = ex.Trials(i).period(u).FREQ;
            ex.Trials(i).period(u).(['lfp_' bands{b} '_pow']) = ...
                mean(ex.Trials(i).period(u).POW(freq >= range{b}(1) & freq <= range{b}(2)));
        end
    end
end
