function stlfp = getSTA(lfptrace, lfptime, spktime, wnd)
% generic function to get spike-triggered averaging LFP

nspk = length(spktime);
ncol = length(-wnd:0.001:wnd);
stlfp = nan(nspk, ncol);
for i = 1:nspk
    tspk = find(lfptime <= spktime(i), 1, 'last');
    tstrt = tspk - (wnd*1000);
    tend = tstrt + ncol -1;
    try
        stlfp(i,:) = lfptrace(tstrt:tend);
    catch
        continue
    end
end
stlfp(sum(isnan(stlfp), 2)==1, :) = [];