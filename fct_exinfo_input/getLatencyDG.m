function [lat, pval, psth_out, ntrial] = getLatencyDG(exinfo, oTrials, smooth_flag)
% [lat, pval, psth_out, ntrial] = getLatencyDG(exinfo, oTrials)
%  or
% [lat, pval, psth_out, ntrial] = getLatencyDG(exinfo, oTrials, smooth_flag)
% 
% returns the psth and the thereof computed response latency estimate
% response and its p-value. the latency is computed using ML-estimation on
% the smoothed psth (Friedmann and Priebe, 1998). the psth is computed for
% the stimulus triggered spike trains in oTrials.
%
% If smooth_flag is true, psth_out is the smoothed version.
% 
% @CL


ntrial = length(oTrials); % number of trials
kernel = ones(50,1)/50; % 50ms boxcar used to smooth the psth

psth = getPSTH(oTrials); % the PSTH
psth_smooth = filter(kernel, 1, psth); % the smoothed PSTH


if ~isempty(oTrials)
    sprintf('id %d par %1.1f', exinfo.id, oTrials(1).(exinfo.param1));
    [lat, ~, pval] = friedmanpriebe(psth_smooth, 'minTheta', 20, ...
        'maxTheta', 100, 'responseSign', 1);
else
    lat = 0;
    pval = 1;
end


% refine the output
psth_out = psth;
if nargin == 3 
    if smooth_flag 
        psth_out = psth_smooth;
    end
end

end