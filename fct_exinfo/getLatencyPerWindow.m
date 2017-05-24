function [latmx, psth] = getLatencyPerWindow(exinfo, ex)
% returns a matrix with response latencies for each stimulus condition
% divided for trial that came first in the fixation period and trials that
% came later (2sec trials with 4 stimuli). psth is an additional cell
% containing the PSTHs the corresponding number of trials.
% 
% @CL

% extent the Trial structure
oTrials = ex.Trials([ex.Trials.Reward]==1); % only valid trials
oTrials = addWindow(exinfo, oTrials); % add window number

% if it is OR, collapse stimuli orientation
if strcmp(exinfo.param1, 'or')
    ind = [oTrials.or]<1000;
    stimor = mod([oTrials(ind).or], 180); stimor = num2cell(stimor);
    [oTrials(ind).or] = deal(stimor{:});
end

% stimuli without blank
stimvls = unique( [ oTrials.(exinfo.param1) ] );
stimvls = stimvls( stimvls < 1000 );


% define the output variables
latmx = nan(length(stimvls), 4);
psth = cell(length(stimvls), 4);


% all trials that are the first within the fixation period
for stim_i = 1:length(stimvls)
    
    ind = [oTrials.(exinfo.param1)]==stimvls(stim_i) & [oTrials.window]==1;
    
    [latmx(stim_i, 1), psth{stim_i, 1, 1}] = getLatencyDG( exinfo, oTrials(ind));
    psth{stim_i, 1, 2} = sum(ind);
    
end


% others - all following trials can be grouped
for stim_i = 1:length(stimvls)
    
    ind = [oTrials.(exinfo.param1)]==stimvls(stim_i) & [oTrials.window]>1;
    
    [latmx(stim_i, 2), psth{stim_i, 2, 1}] = getLatencyDG( exinfo, oTrials(ind));
    psth{stim_i, 2, 2} = sum(ind);
    
end



end
