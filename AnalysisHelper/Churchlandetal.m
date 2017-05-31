function [ spikes ] = Churchlandetal( Trials )
%Churchlandetal
% Trials without reward are excluded. Each neuron in each condition is one
% date point.
%

addTime = 0.200;

for i = 1:length(Trials)
    
    % frame presentation times
    t_frame     = Trials(i).Start - Trials(i).TrialStart;
    if Trials(i).adapt
        t_frame = t_frame(t_frame > Trials(i).adaptationDuration);
    end
    
    t_start     = t_frame(1);
    t_end       = t_frame(end) + mean(diff(t_frame));
    t_before    = t_start - addTime;
    t_after     = t_end + addTime;
    
    
    spk_before = round((Trials(i).Spikes(Trials(i).Spikes >= t_before &...
        Trials(i).Spikes <= t_start) -t_before)  *1000)+1;
    
    spk_pres = round((Trials(i).Spikes(Trials(i).Spikes > t_start &...
        Trials(i).Spikes <= t_end) -t_before) *1000)+1;
    
    spk_after = round((Trials(i).Spikes(Trials(i).Spikes > t_end &...
        Trials(i).Spikes <= t_after) -t_before) *1000)+1;
    
    
    % if the trial length is shorter than 460 ms
    % the time after trial is postponed, e.g. addtional
    % entries with 0 are added in the time between
    dur_pres = round((t_end-t_start) *1000);
    if dur_pres < 460
        diff_ = 460 - dur_pres;
        spk = [spk_before', spk_pres', spk_after'+diff_];
    elseif dur_pres > 460
        long_flag = 1;
        diff_ = dur_pres - 460;
        spk = [spk_before', spk_pres(spk_pres <= 460)', spk_after'-diff_];
    else
        spk = [spk_before', spk_pres', spk_after'];
    end
     
    
    spikes(i, 1:(460+1+2000*addTime)) = false;
    spikes(i, spk) = true;
    
end



end