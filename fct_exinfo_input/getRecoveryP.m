function exinfo = getRecoveryP( exinfo )
% add a field that points out the probability to have a recovered control

clc
% get the mean bootstrapped data
for i = 1:length(exinfo)
    if ~exinfo(i).isRC;     
        exinfo(i) = bootstrap_exinfo( exinfo(i), true ); 
    end
        
        exinfo(i).ret2base = false;
        exinfo(i).iscmp2base = true;
        exinfo(i).fig_recovery = [];
end



% loop through each unit and get the bootstrapped resampling and the
% results whether it is classified as returned to baseline
uq_id = unique([exinfo.id]);

for id = uq_id
    
    if min([exinfo([exinfo.id]==id).dt_id]) >0
        fprintf('this is not the truely first experiment %1.1f\n', exinfo(find([exinfo.id]==id, 1, 'first')).id);
    end
    
    idx = [exinfo.id]==id & ~[exinfo.isRC];
    exinfo(idx) = getRecID(exinfo(idx), true);
    
end

end

function exinfo = getRecID(exinfo, clustercheck)
% groups the experiments according to stimulus and cluster
if clustercheck 
  
    exinfo([exinfo.isc2]) = getRecID(exinfo([exinfo.isc2]), false);
    exinfo(~[exinfo.isc2]) = getRecID(exinfo(~[exinfo.isc2]), false);

elseif ~isempty(exinfo)
    uq_stim = unique({exinfo.param1});
    
    for stim = uq_stim
       idx = strcmp({exinfo.param1}, stim);
       exinfo(idx) = getRecStim(exinfo(idx), true); 
    end
end
end


function exinfo = getRecStim(exinfo, checkocul)
% compares the baseline and control bootstrapped mean tc difference to 0

if checkocul
    exinfo([exinfo.ocul]==1) = getRecStim(exinfo([exinfo.ocul]==1), false);
    exinfo([exinfo.ocul]==0) = getRecStim(exinfo([exinfo.ocul]==0), false);
    exinfo([exinfo.ocul]==-1) = getRecStim(exinfo([exinfo.ocul]==-1), false);
    exinfo([exinfo.ocul]==-11) = getRecStim(exinfo([exinfo.ocul]==-11), false);
    
elseif ~isempty(exinfo)
    
    if exinfo(1).dt_id >0
        iscmp2base = false;
    else
        iscmp2base = true;
    end
    
    [~, firstidx] = min([exinfo.dt_id]);
    mu = exinfo(firstidx).nrep' * exinfo(firstidx).ratemn / sum(exinfo(firstidx).nrep);
    
    
    for k = 1:length(exinfo)
    
        ct = prctile(exinfo(k).rate_resmpl, [5 95]);
        exinfo(k).ret2base = mu > ct(1) && mu < ct(2);
                
        exinfo(k).iscmp2base = iscmp2base;
        

        
        
    end
%     exinfo = plotRecoveryAndEffect(exinfo);
end

end



function exinfo = plotRecoveryAndEffect(exinfo)
% plot the 90% confidence interval for each experiment

fname = [exinfo(1).figname(1:18) exinfo(1).param1 getOcul(exinfo(1).ocul)];
h = figure('Name', fname);


% mean firing rate and confidence intervals
axes('Position', [0.1 0.1 0.8 0.6]);
x = 1:length(exinfo)*2;
fill([x fliplr(x)], [repmat(exinfo(1).mntc_CI_base(1), length(x)) ...
    repmat(exinfo(1).mntc_CI_base(end), length(x))], ...
    'k', 'FaceAlpha', '0.05', 'EdgeAlpha', '0.05');
hold on;

j = 1;
k = 1;

while j <= length(x)
    col = getCol(exinfo(k));
    
    wmu = exinfo(k).ratemn' * exinfo(k).nrep / sum(exinfo(k).nrep);
    plot([j j], exinfo(k).mntc_CI_base([1,5]), 'k'); hold on;
    scatter(j, wmu, col, 'ButtonDownFcn', { @OpenExPlot, exinfo(k) } );
    
    wmu_drug = exinfo(k).ratemn_drug' * exinfo(k).nrep_drug / sum(exinfo(k).nrep_drug);
    plot([j j]+1, exinfo(k).mntc_CI_drug([1,5]), 'k'); hold on;
    scatter(j+1, wmu_drug, col, 'filled', 'ButtonDownFcn', { @OpenExPlot, exinfo(k) } );
    
    k = k+1;
    j = j+2;
end

xlabel('experiment number'); ylabel('averaged TC +/- 5% percentiles');
xlim([.9 max(x)+0.1]);


% information regarding recovery, dose and time of recording
axes('Position', [0.1 0.7 0.8 0.1]);

k = 1; j = 1;
while j <= length(x)
    text(j, 1, sprintf('recovp = %1.2f \ndose = %1.0fnA \nt2strt = %1.2f \n me=%1.0f', ...
        exinfo(k).ret2base, exinfo(k).dose, exinfo(k).dt_id, exinfo(k).ocul), 'FontSize', 8);
    k = k+1;
    j = j+2;
end

xlim([.9 max(x)+0.1]);
axis off

fname = ['C:\Users\Corinna\Documents\Code\plots\Recovery\' fname];
set(h, 'UserData', exinfo);
savefig(h, fname);
delete(h);


for k = 1:length(exinfo)
    exinfo(k).fig_recovery = fname;
end


end



function OpenExPlot(~, ~, exinfo)

openfig(exinfo.fig_tc);
end


function ocul = getOcul(me)

switch me
    case -1
        ocul = '_left';
    case 0
        ocul = '_both';
    case 1
        ocul = '_right';
    case -11
        ocul = '_unknown';
end
end