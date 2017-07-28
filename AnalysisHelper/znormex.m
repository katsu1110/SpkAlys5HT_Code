function [ex, spkrate, spkcount] = znormex(ex, exinfo)
%[ex, spkrate, spkcount] = znormex(ex, exinfo)
% 
% computes the stimulus elicited spike count and resulting spike rate and
% saves the z-normed activity in ex.Trials  for each stimulus.
% 
% additionally, cell arrays with information for the raster plot are saved
% in ex (ex.raster, ex.rastercum, ex.raster, ex.trials_n). The cells are
% sorted corresponding to spkrate.x with x being the stimulus feature (or,
% co, sz, sf,...)
% 
% @CL


param1 = exinfo.param1;
parvls = unique( [ ex.Trials.(param1) ] );

j = 1;
Fs = 1000; % recording frequency

if exinfo.isadapt
    time = 0:1/Fs:5;
else
    time = 0.001:1/Fs:0.45;
end


for par = parvls
    
    % z-scored spikes for this stimulus
    ind = par == [ ex.Trials.(param1) ];
    spkrate(j).nrep = sum(ind);
    spkrate(j).(param1) = par;

    
    % spike rates statistics
    spkrate(j).mn  = mean([ex.Trials(ind).spkRate]);
    spkrate(j).var = var([ex.Trials(ind).spkRate]);
    spkrate(j).sd = std([ex.Trials(ind).spkRate]);
    spkrate(j).sem = spkrate(j).sd / sqrt( spkrate(j).nrep );
    spkrate(j).raw = [ex.Trials(ind).spkRate];
    spkrate(j).resamples = resample2([ex.Trials(ind).spkRate]);
  
    % spike count statistics
    spkcount(j).(param1) = par;
    spkcount(j).mn  = mean([ex.Trials(ind).spkCount]);
    spkcount(j).var = var([ex.Trials(ind).spkCount]);
    
    
    % z normed spike counts
    z = zscore( [ ex.Trials(ind).spkCount ] );
    z = num2cell(z);
    [ex.Trials(ind).zspkcount] = deal(z{:});
    
    % raster - contains 0 for times without spike and 1 for times of spikes
    % is reduced for times between -0.5s and 0.5s around stimulus onset
    idxct = find(ind);
    ct = ex.Trials(idxct);    % current trials
    
    ex.trials_n(j) = sum(ind);
    ex.raster{j,1} = zeros(length(ct), length(time));
    ex.rastercum{j,1} = nan(length(ct), length(time));
        
    % convert the spike times into a matrix raster with bits
    for k = 1:length(idxct)
        
        t_strt = ct(k).Start - ct(k).TrialStart;
        if exinfo.isadapt
            t_strt = t_strt(t_strt > ex.stim.vals.adaptationDur);
        end
        
        idx = round((ct(k).Spikes(...
            ct(k).Spikes>=t_strt(1) & ct(k).Spikes<=t_strt(end))-t_strt(1)) ...
            *Fs);
        
        idx(idx==0) = 1; % avoid bad indexing
        
        ex.raster{j}(k, idx) = 1;
        ex.rastercum{j}(k, idx) = k;
        
    end
    
    j = j+1;
    
end


ex.raster_pval = parvls;

end



function res = resample2(A)
% resample 1000 times from A

n = length(A); 
nsmpl = 1000; % number 

if n > 0
res = ones(n, nsmpl); %initialize variable to prevent overhead
idx = randi(n, 1000, n); % combination of indices corresponding to A's data

for i = 1:nsmpl
    res(:, i) = A(idx(i, :)); %assigning the randomized resamples
end

else
   res = []; 
end

end








