function getNExp(exinfo)
clc


exinfo_RC = exinfo([exinfo.isRC]);
exinfo_DG = exinfo(~[exinfo.isRC]);
exinfo_or = exinfo(strcmp({exinfo.param1}, 'or') & ~[exinfo.isRC]);
exinfo_co = exinfo(strcmp({exinfo.param1}, 'co'));
exinfo_sf = exinfo(strcmp({exinfo.param1}, 'sf'));
exinfo_sz = exinfo(strcmp({exinfo.param1}, 'sz'));


minrep = 4;
minspk = 10;
is5HT = [exinfo.is5HT];

writeNExp(exinfo, 'All', minrep, minspk);
writeNExp(exinfo(is5HT), 'All - 5HT', minrep, minspk);
writeNExp(exinfo(~is5HT), 'All - NaCl', minrep, minspk);

% writeNExp(exinfo_DG, 'All DG', minrep, minspk);
% writeNExp(exinfo_DG([exinfo_DG.is5HT]), 'All DG - 5HT', minrep, minspk);
% writeNExp(exinfo_DG(~[exinfo_DG.is5HT]), 'All DG - NaCl', minrep, minspk);
% 
% writeNExp(exinfo_or([exinfo_or.is5HT]), 'OR - 5HT', minrep, minspk);
% writeNExp(exinfo_or(~[exinfo_or.is5HT]), 'OR - NaCl', minrep, minspk);
% 
% writeNExp(exinfo_sf([exinfo_sf.is5HT]), 'SF - 5HT', minrep, minspk);
% writeNExp(exinfo_sf(~[exinfo_sf.is5HT]), 'SF - NaCl', minrep, minspk);
% 
% writeNExp(exinfo_co([exinfo_co.is5HT]), 'CO - 5HT', minrep, minspk);
% writeNExp(exinfo_co(~[exinfo_co.is5HT]), 'CO - NaCl', minrep, minspk);
% 
% writeNExp(exinfo_sz([exinfo_sz.is5HT]), 'SZ - 5HT', minrep, minspk);
% writeNExp(exinfo_sz(~[exinfo_sz.is5HT]), 'SZ - NaCl', minrep, minspk);
% 
% writeNExp(exinfo_RC([exinfo_RC.is5HT]), 'RevCorr - 5HT', minrep, 0.1);
% writeNExp(exinfo_RC(~[exinfo_RC.is5HT]), 'RevCorr - NaCl', minrep, 0.1);



end


function writeNExp(exinfo, exname, repmin, spkmin)

isbr = [exinfo.electrodebroken];
isbrincl = ~[exinfo.electrodebroken_excl];

fprintf(['\n\n------------------' exname '\n']);

fprintf('                                                 all (ma/ka)       broken     broken&incl \n');

% all data
fprintf('#units                                           ');
writeLine(exinfo);
writeLine(exinfo(isbr));
writeLine(exinfo(isbrincl));
writeMethodStats(exinfo(isbrincl));

% 
% % general inclusion criteria
% fprintf('\n#units ANOVA p<1%%, nrep>=4, min 10 spk/s         \n');
% crit = cellfun(@min, {exinfo.nrep})>=repmin &  cellfun(@min, {exinfo.nrep_drug})>=repmin & ...
%      cellfun(@max, {exinfo.ratemn})>=spkmin & ...
%      [exinfo.p_anova]<0.01 & [exinfo.p_anova_drug]<0.01 & ...
%      cellfun(@(x) isnan(x) || x<150, {exinfo.resistance});
% writeLine(exinfo(crit));
% writeLine(exinfo(isbr&crit));
% writeLine(exinfo(isbrincl&crit));

% % regression 
% fprintf('\n#units gen crit and regression R2 > 70%%          ');
% crit2 = crit & [exinfo.r2reg] > 0.7;
% writeLine(exinfo(crit2));
% writeLine(exinfo(isbr&crit2));
% writeLine(exinfo(isbrincl&crit2));
% 
% 
% % tuning curve
% fprintf('\n#units gen crit and tuning R2 > 70%%              ');
% crit3 = crit & [exinfo.gaussr2] > 0.7 & [exinfo.gaussr2_drug] > 0.7;
% writeLine(exinfo(crit3));
% writeLine(exinfo(isbr&crit3));
% writeLine(exinfo(isbrincl&crit3));
% 
% 
% % tuning curve & regression
% fprintf('\n#units gen crit and both fits R2 > 70%%           ');
% crit3 = crit2 & [exinfo.gaussr2] > 0.7 & [exinfo.gaussr2_drug] > 0.7;
% writeLine(exinfo(crit3));
% writeLine(exinfo(isbr&crit3));
% writeLine(exinfo(isbrincl&crit3));


end


function writeLine(exinfo)

ic2 = [exinfo.isc2];
ma = [exinfo.ismango];

n_all = length(unique([exinfo(ic2).id])) + length(unique([exinfo(~ic2).id]));
n_ma = length(unique([exinfo(ic2&ma).id])) + length(unique([exinfo(~ic2&ma).id]));
n_ka = length(unique([exinfo(ic2&~ma).id])) + length(unique([exinfo(~ic2&~ma).id]));
fprintf('%1.0f (%1.0f/%1.0f)        ', n_all, n_ka, n_ma);



% unique([exinfo(ic2).id])
end


function writeMethodStats(exinfo)

fprintf('\n');

%----- variables that are only counted once per unit
unique_id = unique([exinfo.id]);
ids = []; 
for i = 1:length(unique_id)
   ids = [ids, find([exinfo.id] == unique_id(i), 1, 'first')];
end

% inclusion indices
idx = zeros(1, length(exinfo));
idx(ids) = true;
ic2 = [exinfo.isc2];

rf_sz = prctile([exinfo(~ic2 & idx).RFwy], [0,50,100]); % receptive field size (equivalent width)
ecc = prctile([exinfo(~ic2 & idx).ecc], [0,50,100]);  % eccentricity


%----- variables that are merged across units and clusters
dt_drug = prctile([exinfo.dose], [0,50,100]); 

idx_recov = ~ic2 & [exinfo.recov_p]>0.05 & [exinfo.t_recov]>=0;
rec_t = prctile([exinfo(idx_recov).t_recov], ...
    [0,50,100]);  % duration of recovery for fully recovered data 


fprintf(['median and range for all broken&incl units... \n' ...
    ' RF size (=eqw) (1 per unit):          %1.2f [%1.2f/%1.2f] \n'...
    ' eccentricity (1 per unit):            %1.2f [%1.2f/%1.2f] \n'...
    ' recovery time [s] (fully rec. only):  %1.2f [%1.2f/%1.2f] \n'...
    '   n recovered:                        %1d'], ...
    rf_sz([2,1,3]), ecc([2,1,3]), rec_t([2,1,3]), sum(idx_recov));



end