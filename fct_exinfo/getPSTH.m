function [psth, t] = getPSTH(Trials)
% returns the psth normalized to spks/s and corresponding time vector
% 
% takes a struct of trials from a ex.mat file and converts the
% spike trains to a psth. Time bin width is 1ms.
%
% @Corinna Lorenz, 24.08.2016

fs = 1000; % sampling frequency

allspk = []; 
for i = 1:length(Trials)
    % preallocate variables
    frameon = Trials(i).Start - Trials(i).TrialStart;   % frame start times
    ts(i) = frameon(1);                                 % first frame                      
    te(i) = frameon(end) + mean(diff(frameon));         % last frame ending

    % spikes within the presentation time
    spks = Trials(i).Spikes ( Trials(i).Spikes >= ts(i) & Trials(i).Spikes <= te(i) );
    allspk = [allspk; spks-ts(i)] ;    % align to stimulus onset
end

% time bin vector
t = 1/fs: 1/fs :round(mean(te-ts)*1000)/1000;

% compute psth
psth = histc(allspk, t); 
psth = psth ./ length(Trials) .*fs; %convert to spk/s
if isempty(psth); psth=zeros(length(t),1); end
if size(psth,2)>1; psth=psth'; end
if mod(length(psth), 2); psth = psth(1:end-1); t = t(1:end-1); end

