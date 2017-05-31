function exinfo = getPhaseSelectivityRC( exinfo )




for i = 1:length(exinfo)
    
    % only use flashed grating experiments
    if ~exinfo(i).isRC %|| exinfo(i).id ~= 277
        continue;
    end
  

    %load ex file and compute reverse correlation subspace analysis for
    %two dimensions: or and phase
    ex = loadCluster(exinfo(i).fname);
    res = HN_computeLatencyAndNetSpk_PhaseSelectivity([],ex);
    exinfo(i).phasesel = r0r1(res);

    % repeat the commands for the 5HT and NaCl blocks
    ex = loadCluster(exinfo(i).fname_drug);
    res_drug = HN_computeLatencyAndNetSpk_PhaseSelectivity([],ex);
    exinfo(i).phasesel_drug = r0r1(res_drug);
    
    
    
    % plot the results as colored blocks
    figure('name', exinfo(i).figname);
    s(1)= subplot(1,2,1);
    [y, idx] = sort(res.sdfs.y,2);
    imagesc( y(1,:), res.sdfs.x(:, 1), res.netSpikesPerFrame(:, idx(1,:)));
    title(sprintf('Baseline f1/f0 = %1.2f', exinfo(i).phasesel));
    colorbar('southoutside'); xlabel('phase'); ylabel('orientation')
    
    s(2) = subplot(1,2,2);
    [y, idx] = sort(res_drug.sdfs.y,2);
    imagesc( y(1,:), res_drug.sdfs.x(:, 1), res_drug.netSpikesPerFrame(:, idx(1,:)));
    title(sprintf('%s f1/f0 = %1.2f', exinfo(i).drugname, exinfo(i).phasesel_drug));
    colorbar('southoutside');xlabel('phase'); ylabel('orientation')
    
    clim_ = [s(1).CLim, s(2).CLim];
%     set(s, 'CLim', [min(clim_), max(clim_)]);
    savefig(gcf, exinfo(i).fig_phase);
    delete(gcf);
end


end



function phasesel = r0r1(res)
% compute the phase selectivity, i.e. f0/f1 

% get the components from frames with preferred orientation only
[~, maxi] = max( mean( res.netSpikesPerFrame, 2 ) ); 
spks_prefor = res.netSpikesPerFrame(maxi, :);

% negative spike rates falsify the phase selectivity. offset all to avoid the problem.
% if any(spks_prefor<0)
%     spks_prefor = spks_prefor + abs(min(spks_prefor));
% end

% f0 is the mean firing rate
f0 = mean(spks_prefor) *2;

% f1 is the amplitude of the first harmonic
% I compute it as the difference between maximum and minimum phase response
f1 = max(spks_prefor)-min(spks_prefor) /2;

phasesel = f1/f0;

end