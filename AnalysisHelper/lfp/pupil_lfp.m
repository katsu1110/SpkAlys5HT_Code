function [pslfp] = pupil_lfp(exinfo)
% examine correlation between pupil size and lfp power, with or without 5HT
% INPUT: exinfo ... single pair of experiment
% OUTPUT: pslfp ... output structure

%%
% structure initialization
pslfp.pupil_timecourse = [];
pslfp.drugname = exinfo.drugname;

% load data ---------------------------------
% control
ex0 = loadCluster(exinfo.fname, 'ocul', exinfo.ocul, 'loadlfp', true);
% drug
ex2 = loadCluster(exinfo.fname_drug, 'ocul', exinfo.ocul, 'loadlfp', true);

% further preprocessing of LFP -----------------
[ex0] = filterLFP(ex0);
[ex2] = filterLFP(ex2);

% spike-triggered LFP --------------------------------
[~, ~, ~, ~, ~, ~, ex0] = spktriglfp(ex0);
[~, ~, ~, ~, ~, ~, ex2] = spktriglfp(ex2);

% preprocess pupil data ----------------------------
[ex0_sps, ex0_lps, ex0] = pupilSplit(ex0);
[ex2_sps, ex2_lps, ex2] = pupilSplit(ex2);

% reverse correlation analysis
if exinfo.isRC
    for i = 1:4
        switch i
            case 1
                exd = ex2_sps;
                labd = ['pupil small, ' pslfp.drugname];
            case 2
                exd = ex2_lps;
                labd = ['pupil large, ' pslfp.drugname];
            case 3
                exd = ex0_sps;
                labd = 'pupil small, baseline';
            case 4
                exd = ex0_lps;
                labd = 'pupil large, baseline';            
        end               
        [stmMat, actMat] = ex4RCsub(exd, 'or', 'LFP_prepro');
        pslfp.rcsub(i).results = reverse_corr_subspace(stmMat, actMat, 300, 100, 0);
        pslfp.rcsub(i).label = labd;
        pslfp.inter_table_lat(i) = pslfp.rcsub(i).results.latency;
    end
    % gain or additive change (normalized by large ps max)
    m = max([pslfp.rcsub(4).results.stm.peak]);
    pslfp.type2reg.base = gmregress([pslfp.rcsub(3).results.stm.peak]/m, [pslfp.rcsub(1).results.stm.peak]/m);
    m = max([pslfp.rcsub(2).results.stm.peak]);
    pslfp.type2reg.drug = gmregress([pslfp.rcsub(4).results.stm.peak]/m, [pslfp.rcsub(2).results.stm.peak]/m);
end


% trial number
len_tr0 = length(ex0.Trials);
len_tr2 = length(ex2.Trials);    

% store pupil size time-course ---------------------
% -- drug -- label_tr -- TC ---
nmin = 1000;
for i = 1:len_tr0
    if length(ex0.Trials(i).pupil_z) < nmin
        nmin = length(ex0.Trials(i).pupil_z);
    end
end
for i = 1:len_tr2
    if length(ex2.Trials(i).pupil_z) < nmin
        nmin = length(ex2.Trials(i).pupil_z);
    end
end

psmat0 = zeros(len_tr0 , 2 + nmin);
for i = 1:len_tr0
    psmat0(i,2) = ex0.Trials(i).n_stm;
    psmat0(i,3:end) = ex0.Trials(i).pupil_z(1:nmin);
end
psmat2 = ones(len_tr2 , 2 + nmin);
for i = 1:len_tr2
    psmat2(i,2) = ex2.Trials(i).n_stm;
    psmat2(i,3:end) = ex2.Trials(i).pupil_z(1:nmin);
end

pslfp.pupil_timecourse = [psmat0; psmat2];

% use only the last stimulus in the DG experiments
if ~exinfo.isRC
    ex0.Trials = ex0.Trials(psmat0(:,2)==4);
    ex2.Trials = ex2.Trials(psmat2(:,2)==4);
    len_tr0 = sum(psmat0(:,2)==4);
    len_tr2 = sum(psmat2(:,2)==4);
end

% make matrices for GLM ----------------------------------
l = 14;
mat0 = zeros(len_tr0, l);
mat2 = ones(len_tr2, l);

% baseline
for i = 1:len_tr0        
    % lfp powers
    mat0(i,2) = ex0.Trials(i).lfp_delta_pow;
    mat0(i,3) = ex0.Trials(i).lfp_theta_pow;
    mat0(i,4) = ex0.Trials(i).lfp_alpha_pow;
    mat0(i,5) = ex0.Trials(i).lfp_beta_pow;
    mat0(i,6) = ex0.Trials(i).lfp_gamma_pow;
    
    % stlfp powers
    mat0(i,7) = ex0.Trials(i).stlfp_delta;
    mat0(i,8) = ex0.Trials(i).stlfp_theta;
    mat0(i,9) = ex0.Trials(i).stlfp_alpha;
    mat0(i,10) = ex0.Trials(i).stlfp_beta;
    mat0(i,11) = ex0.Trials(i).stlfp_gamma;
    
    % average lfp
    mat0(i,12) = nanmean(ex0.Trials(i).LFP_prepro(ex0.Trials(i).LFP_prepro_time > 0.35));
    
    % stLFP amplitude (t=0)
    mat0(i, 13) = nanmin(ex0.Trials(i).mean_stLFP);

    % pupil size
    mat0(i, l) = ex0.Trials(i).pupil_val;
end

% drug
for i = 1:len_tr2        
    % lfp powers
    mat2(i,2) = ex2.Trials(i).lfp_delta_pow;
    mat2(i,3) = ex2.Trials(i).lfp_theta_pow;
    mat2(i,4) = ex2.Trials(i).lfp_alpha_pow;
    mat2(i,5) = ex2.Trials(i).lfp_beta_pow;
    mat2(i,6) = ex2.Trials(i).lfp_gamma_pow;

    % stlfp powers
    mat2(i,7) = ex2.Trials(i).stlfp_delta;
    mat2(i,8) = ex2.Trials(i).stlfp_theta;
    mat2(i,9) = ex2.Trials(i).stlfp_alpha;
    mat2(i,10) = ex2.Trials(i).stlfp_beta;
    mat2(i,11) = ex2.Trials(i).stlfp_gamma;
    
    % average lfp
    mat2(i,12) = nanmean(ex2.Trials(i).LFP_prepro(ex2.Trials(i).LFP_prepro_time > 0.35));
    
    % stLFP amplitude (presumably t=0)
    mat2(i, 13) = nanmin(ex2.Trials(i).mean_stLFP);
    
    % pupil size
    mat2(i, l) = ex2.Trials(i).pupil_val;
end

% correlation ------------
for d = 1:2
    switch d
        case 1
            mat = mat0;
            fieldname = 'control';
        case 2
            mat = mat2;
            fieldname = 'drug';
    end
    for i = 1:l-2
        [pslfp.corr.(fieldname).rho(i), pslfp.corr.(fieldname).pval(i)] = corr(mat(:,i+1), mat(:, l), 'type', 'Spearman');
    end
end

% interaction table --------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%
% LFP %% 5HT %%% base %%%

% S-ps      %%%      %%%
% L-ps      %%%      %%%
%%%%%%%%%%%%%%%%%%%%%%%%

med0 = median(mat0(:, l));
med2 = median(mat2(:, l));
for i = 1:l-2
    pslfp.interaction(i).table(1,1) = nanmean(mat2(mat2(:, l) < med2, i+1));
    pslfp.interaction(i).table(2,1) = nanmean(mat2(mat2(:, l) > med2, i+1));
    pslfp.interaction(i).table(1,2) = nanmean(mat0(mat0(:, l) < med0, i+1));
    pslfp.interaction(i).table(2,2) = nanmean(mat0(mat0(:, l) > med0, i+1));
    
    % normalize the table by the grand mean
    pslfp.interaction(i).table = pslfp.interaction(i).table/...
        mean(mean(pslfp.interaction(i).table));
end
