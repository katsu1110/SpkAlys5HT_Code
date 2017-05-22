function [lat, pval, psth_out, ntrial] = getLatencyDG(exinfo, oTrials, smooth_flag)
% returns latency and unsmoothed, but normalized psth for trials with
% drifting grating



ntrial = length(oTrials);
psth = getPSTH(oTrials);


% psth smoothing
kernel = ones(50,1)/50;

psth_smooth = filter(kernel, 1, psth);


if ~isempty(oTrials)
    sprintf('id %d par %1.1f', exinfo.id, oTrials(1).(exinfo.param1));
    [lat, ~, pval] = friedmanpriebe(psth_smooth, 'minTheta', 20, ...
        'maxTheta', 100, 'responseSign', 1);
else
    lat = 0;
    pval = 1;
end




% output
psth_out = psth;
if nargin == 3 
    if smooth_flag 
        psth_out = psth_smooth;
    end
end

end