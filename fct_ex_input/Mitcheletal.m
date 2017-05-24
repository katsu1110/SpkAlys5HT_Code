function [ mitchel_mn, mitchel_var] = Mitcheletal( Trials, param1 )
%Mitchelletal
% Trials without reward are excluded. Each neuron in each condition is one
% date point.
%


nbin = 4;

Trials  = getSpkPerBin(Trials, nbin);
mSpk    = getMeanSpkPerCond(Trials, param1);

mitchel_mn = horzcat(mSpk.mn);
mitchel_var= horzcat(mSpk.var);

end



function Trials = getSpkPerBin(Trials, nbin)
% returns spikes binned for every 100ms for each trial

for n = 1:length(Trials)
    
    t_spks = Trials(n).Spikes; % time of spike occurance
    trial_strt = Trials(n).Start - Trials(n).TrialStart; % frame start times
    if Trials(n).adapt
        trial_strt = trial_strt(trial_strt > Trials(n).adaptationDuration);
    end
    
    trial_strt(end+1) = trial_strt(end) + mean(diff(trial_strt)); % time of last stimulus frame offset
    
    
    tstep = 0.1; % steps of 100ms
    offset = 0.05; % start binning after 50ms after stimulus onset
    for bin_i = 1:nbin
        bin_strt(bin_i)    = tstep*(bin_i-1)+trial_strt(1)+offset; 
        bin_end(bin_i)     = bin_strt(bin_i)+tstep;
        
        if bin_end(bin_i) > trial_strt(end)
            bin_end(bin_i) = trial_strt(end); 
        end
        spkin100msbin(bin_i) = sum( t_spks>=bin_strt(bin_i) & t_spks<bin_end(bin_i) ) ;
        
    end
        
    Trials(n).spkin100msbin = spkin100msbin;
end

end


function mSpk = getMeanSpkPerCond(Trials, param)
% evaluates mean and variance of each time bin over all neurons and
% % conditions
% from Mitchell et al. (2007), analysis
% The Fano factor was computed in nonoverlapping 100 ms bins. The spike
% count was computed in each bin for each trial. Then, for each bin, the
% Fano factor was computed as the ratio of variance in spike counts across
% trials to spike count mean.


% stimulus parameter
parVal = unique([Trials.(param)]);
    
for i = 1:length(parVal)
    
    idx = [Trials.(param)] == parVal(i);
    mSpk(i).mn = mean(vertcat(Trials(idx).spkin100msbin), 1);
    mSpk(i).var = var(vertcat(Trials(idx).spkin100msbin), 0, 1);
    mSpk(i).parVal = parVal(i);
    
end

end

