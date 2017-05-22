function getNExp(exinfo)



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

writeNExp(exinfo_DG, 'All DG', minrep, minspk);
writeNExp(exinfo_DG([exinfo_DG.is5HT]), 'All DG - 5HT', minrep, minspk);
writeNExp(exinfo_DG(~[exinfo_DG.is5HT]), 'All DG - NaCl', minrep, minspk);

writeNExp(exinfo_or([exinfo_or.is5HT]), 'OR - 5HT', minrep, minspk);
writeNExp(exinfo_or(~[exinfo_or.is5HT]), 'OR - NaCl', minrep, minspk);

writeNExp(exinfo_sf([exinfo_sf.is5HT]), 'SF - 5HT', minrep, minspk);
writeNExp(exinfo_sf(~[exinfo_sf.is5HT]), 'SF - NaCl', minrep, minspk);

writeNExp(exinfo_co([exinfo_co.is5HT]), 'CO - 5HT', minrep, minspk);
writeNExp(exinfo_co(~[exinfo_co.is5HT]), 'CO - NaCl', minrep, minspk);

writeNExp(exinfo_sz([exinfo_sz.is5HT]), 'SZ - 5HT', minrep, minspk);
writeNExp(exinfo_sz(~[exinfo_sz.is5HT]), 'SZ - NaCl', minrep, minspk);

writeNExp(exinfo_RC([exinfo_RC.is5HT]), 'RevCorr - 5HT', minrep, 0.1);
writeNExp(exinfo_RC(~[exinfo_RC.is5HT]), 'RevCorr - NaCl', minrep, 0.1);



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


% general inclusion criteria
fprintf('\n#units ANOVA p<1%%, nrep>=4, min 10 spk/s         ');
crit = cellfun(@min, {exinfo.nrep})>=repmin &  cellfun(@min, {exinfo.nrep_drug})>=repmin & ...
     cellfun(@max, {exinfo.ratemn})>=spkmin & ...
     [exinfo.p_anova]<0.01 & [exinfo.p_anova_drug]<0.01 & ...
     cellfun(@(x) isnan(x) || x<150, {exinfo.resistance});
writeLine(exinfo(crit));
writeLine(exinfo(isbr&crit));
writeLine(exinfo(isbrincl&crit));


% regression 
fprintf('\n#units gen crit and regression R2 > 70%%          ');
crit2 = crit & [exinfo.r2reg] > 0.7;
writeLine(exinfo(crit2));
writeLine(exinfo(isbr&crit2));
writeLine(exinfo(isbrincl&crit2));


% tuning curve
fprintf('\n#units gen crit and tuning R2 > 70%%              ');
crit3 = crit & [exinfo.gaussr2] > 0.7 & [exinfo.gaussr2_drug] > 0.7;
writeLine(exinfo(crit3));
writeLine(exinfo(isbr&crit3));
writeLine(exinfo(isbrincl&crit3));


% tuning curve & regression
fprintf('\n#units gen crit and both fits R2 > 70%%           ');
crit3 = crit2 & [exinfo.gaussr2] > 0.7 & [exinfo.gaussr2_drug] > 0.7;
writeLine(exinfo(crit3));
writeLine(exinfo(isbr&crit3));
writeLine(exinfo(isbrincl&crit3));


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