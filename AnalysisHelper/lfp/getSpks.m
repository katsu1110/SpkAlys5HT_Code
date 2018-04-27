function [spk, spkc] = getSpks(trials, wnd)
% spikes within the stimulus presentation time
% (automatically exclude the first 350ms)
% INPUT: trials ... ex.Trials
%              wnd ... [start, end] of analysis window

% time offset for spikes (300ms after stimulus onset as default)
if nargin<2; wnd = [0.35, 0]; end
awnd_strt = wnd(1); %<- 350ms after stimulus onset (default)
awnd_end = wnd(2); %<- 000ms before stimulus end (default)
N = length(trials);
spk = cell(1, N);
spkc = zeros(1, N);
for i = 1:N
    t_strt = trials(i).Start - trials(i).TrialStart;
    t_end = t_strt(end)+mean(diff(t_strt));

    spk_tr = trials(i).Spikes(trials(i).Spikes >= t_strt(1)+awnd_strt & ...
        trials(i).Spikes <= t_end-awnd_end ) - t_strt(1);
    spk{i} = round(spk_tr*1000)/1000;
    spkc(i) = length(spk_tr);
end
