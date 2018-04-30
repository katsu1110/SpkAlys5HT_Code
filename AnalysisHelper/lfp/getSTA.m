function stlfp = getSTA(lfptrace, lfptime, spktime, wnd, fs)
% generic function to get spike-triggered averaging LFP

nspk = length(spktime);
ncol = length(-wnd:(1/fs):wnd);
stlfp = nan(nspk, ncol);
for i = 1:nspk
    [~, tspk] = min(abs(lfptime - spktime(i)));
%     tspk = find(lfptime <= spktime(i), 1, 'last');
    tstrt = tspk - wnd*fs;
    tend = tstrt + ncol - 1;
    try
        stlfp(i,:) = lfptrace(tstrt:tend);
    catch
        continue
    end
end
stlfp(any(isnan(stlfp), 2), :) = [];