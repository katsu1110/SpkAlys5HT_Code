function [latmx, psth] = getLatencyPerWindow(exinfo, ex)
% takes the ex info and fname and returns a matrix with latency for each
% stimuli condition and window of trial sequence (2sec trials with 4
% stimuli). nmx is the additional matrix containin the number of trials in
% the entry.


oTrials = addWindow(exinfo, ex.Trials); % add window number
oTrials = oTrials([oTrials.Reward]==1); % only valid trials

% if it is OR, collapse stimuli orientation
if strcmp(exinfo.param1, 'or')
    ind = [oTrials.or]<1000;
    stimor = mod([oTrials(ind).or], 180); stimor = num2cell(stimor);
    [oTrials(ind).or] = deal(stimor{:});
end


% stimuli parameters without blank
parvls = unique( [ oTrials.(exinfo.param1) ] );
parvls = parvls( parvls < 1000 );


% preset output variables
latmx = nan(length(parvls), 4);
psth = cell(length(parvls), 4);


% first window
for par_i = 1:length(parvls)
    
    ind = [oTrials.(exinfo.param1)]==parvls(par_i) & [oTrials.window]==1;
    
    [latmx(par_i, 1), psth{par_i, 1, 1}] = getLatencyDG( exinfo, oTrials(ind));
    psth{par_i, 1, 2} = sum(ind);
    
end


% others
for par_i = 1:length(parvls)
    
    ind = [oTrials.(exinfo.param1)]==parvls(par_i) & [oTrials.window]>1;
    
    [latmx(par_i, 2), psth{par_i, 2, 1}] = getLatencyDG( exinfo, oTrials(ind));
    psth{par_i, 2, 2} = sum(ind);
    
end



end
