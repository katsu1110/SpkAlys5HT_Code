function [h, list_inc] = plotFigureNCFF(exinfo, nanminss, name, saveoption, varargin)
%% plot scatter for NC and FF for revision
% INPUT: exinfo ... given by Corinna's 'runExinfoAnalysis'
%             nanminss ... 8 or 6 or 4 (the nanminimum sample size per stm)
%             name ... 'nc' or 'ff' or 'ff_peak-blank'
%             
%
% written by Katsuhisa (19.07.17)
% ++++++++++++++++++++++++++++++++++++

close all

switch nargin
    case 0
        error('provide exinfo!')
    case 1
        nanminss = 4;
        name = 'nc';
        saveoption = 0;
    case 2
        name = 'nc';
        saveoption = 0;
    case 3
        saveoption = 0;
end

% basic figure parameters
sz = 10;
a = 0.5;
aa = 0.8;

% set axis range
switch name
    case 'nc'
        xl = [-1 1];
        yl = [-1 1];
        xl2 = [0.06 8];
        yl2 = yl;
        field = 'nc_re';
    case 'ff'
        xl = [0 15];
        yl = [0 15];
        xl2 = [0.06 8];
        yl2 = [-10 10];
        field = 'ff_re';
    case 'ff_peak-blank'
        xl = [-25 25];
        yl = [-25 25];
        xl2 = [0.06 8];
        yl2 = [-10 10];
        field = 'ff_re';
end

% idx for the nanminimum sample size per stm
row = nanminss/2 - 1;


% load list of index
% load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stimulus_cond.mat')
load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stim_cond_2007.mat')
% load('Z:\Corinna\SharedCode\Katsu\incl_i_one_stimulus_cond.mat')
load('Z:\Corinna\SharedCode\Katsu\listpvalue_modulation.mat')

% split exinfo based on the list
for i = 1:length(exinfo)
    exinfo(i).p_mod = p_modulation{i}(2);
    if exinfo(i).param1=='or'
        exinfo(i).stmtype = 1;
    elseif exinfo(i).param1=='sf'
        exinfo(i).stmtype = 2;
    elseif exinfo(i).param1=='co'
        exinfo(i).stmtype = 3;
    elseif exinfo(i).param1=='sz'
        exinfo(i).stmtype = 4;
    end        
end
exinfo = exinfo(incl_i);

% loop to go through all the stimuli types
T = nan(43,15);
base_mango_5ht_all = [];
base_mango_nacl_all = [];
base_kaki_5ht_all = [];
base_kaki_nacl_all = [];
p_mango_5ht_all = [];
p_mango_nacl_all = [];
drug_mango_5ht_all = [];
drug_mango_nacl_all = [];
drug_kaki_5ht_all = [];
drug_kaki_nacl_all = [];
p_kaki_5ht_all = [];
p_kaki_nacl_all = [];
ss_mango_5ht_all = [];
ss_drug_mango_5ht_all = [];
ss_mango_nacl_all = [];
ss_drug_mango_nacl_all = [];
ss_kaki_5ht_all = [];
ss_drug_kaki_5ht_all = [];
ss_kaki_nacl_all = [];
ss_drug_kaki_nacl_all = [];
rr_mango_5ht_all = [];
rr_mango_nacl_all = [];
rr_kaki_5ht_all = [];
rr_kaki_nacl_all = [];
ssperexp = [];
ssperstm = [];
list_inc = [];
h(1) = figure(1);
plot(xl, yl, '-', 'color', 0.6*ones(1,3))
if ~strcmp(name, 'ff')
    hold on;
    plot([0 0],yl,'-', 'color', 0.6*ones(1,3), 'linewidth',0.5)
    hold on;
    plot(xl,[0 0],'-', 'color', 0.6*ones(1,3), 'linewidth',0.5)
end
h(2) = figure(2);
plot([1 1],yl2,'-', 'color', 0.6*ones(1,3), 'linewidth',0.5)
hold on;
plot(xl2,[0 0],'-', 'color', 0.6*ones(1,3), 'linewidth',0.5)

for i = 1:4
    switch i
        case 1
            param = 'or';
            stmtype = 1;
        case 2
            param = 'sf';
            stmtype = 2;
        case 3
            param = 'co';
            stmtype = 3;
        case 4
            param = 'sz';
            stmtype = 4;
    end
    
    disp('++++++++++++++++++++++++++++')
    disp(param)
    
    list_temp = find([exinfo.stmtype]==stmtype);
    exinfo_temp = exinfo(list_temp);
    len_exp = length(exinfo_temp);
    
    % sample size per experiment
    ss = nan(1, len_exp);
    ss_drug = nan(1, len_exp);

    % parameter
    base = nan(1, len_exp);
    drug = nan(1, len_exp);
    
    switch name
        case 'nc'
            disp('storing data of nc...')
            for k = 1:len_exp
                try
                    if isempty(strfind(exinfo_temp(k).fname, 'all.grating')) && ...
                            isempty(strfind(exinfo_temp(k).fname_drug, 'all.grating'))
                        
                        base(k) = exinfo_temp(k).(field)(row, 1);
                        drug(k) = exinfo_temp(k).([field '_drug'])(row, 1);
                        
                        % correction by regression
                        base(k) = base(k) - 0.0062698*[exinfo_temp(k).c0geomn];
                        drug(k) = drug(k) - 0.0062698*[exinfo_temp(k).c0geomn_drug];
                        
                        ss(k) = nansum(exinfo_temp(k).sampleSize_raw(row, 1:end-1),2);
                        ss_drug(k) = nansum(exinfo_temp(k).sampleSize_raw_drug(row, 1:end-1),2);      
                        
                        ssperstm = [ssperstm, ...
                            exinfo_temp(k).sampleSize_raw(row, 1:end-1), ...
                            exinfo_temp(k).sampleSize_raw_drug(row, 1:end-1)];
                    else
                        list_temp(k) = nan;
                        disp('concatenated; skipped...')
                    end                   
                catch
                    list_temp(k) = nan;
                    continue
                end
            end
        case 'ff_peak-blank'
            disp('storing data of ff: peak - blank...')
            for k = 1:len_exp                
                try
                    if isempty(strfind(exinfo_temp(k).fname, 'all.grating')) && ...
                            isempty(strfind(exinfo_temp(k).fname_drug, 'all.grating'))

                        ff_base = exinfo_temp(k).(field)(row, exinfo_temp(k).pfi)...
                            - exinfo_temp(k).(field)(row, end);
                        ff_drug = exinfo_temp(k).([field '_drug'])(row, exinfo_temp(k).pfi_drug)...
                            - exinfo_temp(k).([field '_drug'])(row, end);
                        base(k) = ff_base;
                        drug(k) = ff_drug;

                        ss(k) = nansum(exinfo_temp(k).sampleSize_raw(row, [exinfo_temp(k).pfi end]));
                        ss_drug(k) = nansum(exinfo_temp(k).sampleSize_raw_drug(row, [exinfo_temp(k).pfi_drug end]));

                        ssperstm = [ssperstm, ...
                                exinfo_temp(k).sampleSize_raw(row, [exinfo_temp(k).pfi end]), ...
                                exinfo_temp(k).sampleSize_raw_drug(row, [exinfo_temp(k).pfi_drug end])];
                    else
                        list_temp(k) = nan;
                        disp('concatenated; skipped...')
                    end
                catch
                    list_temp(k) = nan;
                    continue
                end
            end
        case 'ff'            
            disp('storing data of ff...')
            for k = 1:len_exp
                try
                    if isempty(strfind(exinfo_temp(k).fname, 'all.grating')) && ...
                            isempty(strfind(exinfo_temp(k).fname_drug, 'all.grating'))
                        
                        ff_base = exinfo_temp(k).(field)(row, 1:end-1);
                        ff_drug = exinfo_temp(k).([field '_drug'])(row, 1:end-1);
                        idx = intersect(find(~isnan(ff_base)), find(~isnan(ff_drug)));
                        base(k) = mean(ff_base(idx));
                        drug(k) = mean(ff_drug(idx));

                        ss(k) = nansum(exinfo_temp(k).sampleSize_raw(row, idx),2);
                        ss_drug(k) = nansum(exinfo_temp(k).sampleSize_raw_drug(row, idx),2);

                        ssperstm = [ssperstm, ...
                                exinfo_temp(k).sampleSize_raw(row, idx), ...
                                exinfo_temp(k).sampleSize_raw_drug(row, idx)];
                    else
                        list_temp(k) = nan;
                        disp('concatenated; skipped...')
                    end
                catch
                    list_temp(k) = nan;
                    continue
                end
            end
    end
    
    % list of experiments
    list_inc = [list_inc list_temp];
    
    % sample size per experiment
    ssperexp = [ssperexp, ss(k), ss_drug(k)];
    
    % relative rate
     rr = [exinfo_temp.nonparam_ratio];
     
    % mango =======================
    % 5HT
    mango_5ht = find([exinfo_temp.ismango]==1 & [exinfo_temp.is5HT]==1);
    base_mango_5ht = base(mango_5ht);
    drug_mango_5ht = drug(mango_5ht);
    rr_mango_5ht = rr(mango_5ht);
    ss_mango_5ht = ss(mango_5ht);
    ss_drug_mango_5ht = ss_drug(mango_5ht);
    p_mango_5ht = ones(1, length(mango_5ht));
    p_mango_5ht([exinfo_temp(mango_5ht).p_mod] >= 0.05) = 0;
    disp(['Mango, 5HT pairs, significant: ' num2str(length(find(p_mango_5ht==1)))])
    disp(['Mango, 5HT pairs, non-significant: ' num2str(length(find(p_mango_5ht==0)))])
    
    base_mango_5ht_all = [base_mango_5ht_all base_mango_5ht];
    drug_mango_5ht_all = [drug_mango_5ht_all drug_mango_5ht];
    p_mango_5ht_all = [p_mango_5ht_all p_mango_5ht];
    ss_mango_5ht_all = [ss_mango_5ht_all ss_mango_5ht];
    ss_drug_mango_5ht_all = [ss_drug_mango_5ht_all ss_drug_mango_5ht];
    rr_mango_5ht_all = [rr_mango_5ht_all, rr_mango_5ht];
    
    % NaCl
    mango_nacl = find([exinfo_temp.ismango]==1 & [exinfo_temp.is5HT]==0);
    base_mango_nacl = base(mango_nacl);
    drug_mango_nacl = drug(mango_nacl);
    rr_mango_nacl = rr(mango_nacl);
    ss_mango_nacl = ss(mango_nacl);
    ss_drug_mango_nacl = ss_drug(mango_nacl);
    p_mango_nacl = ones(1, length(mango_nacl));
    p_mango_nacl([exinfo_temp(mango_nacl).p_mod] >= 0.05) = 0;
    disp(['Mango, NaCl pairs, significant: ' num2str(length(find(p_mango_nacl==1)))])
    disp(['Mango, NaCl pairs, non-significant: ' num2str(length(find(p_mango_nacl==0)))])
    
    base_mango_nacl_all = [base_mango_nacl_all base_mango_nacl];
    drug_mango_nacl_all = [drug_mango_nacl_all drug_mango_nacl];
    p_mango_nacl_all = [p_mango_nacl_all p_mango_nacl];
    ss_mango_nacl_all = [ss_mango_nacl_all ss_mango_nacl];
    ss_drug_mango_nacl_all = [ss_drug_mango_nacl_all ss_drug_mango_nacl];
    rr_mango_nacl_all = [rr_mango_nacl_all, rr_mango_nacl];
    
    % plot -------------------------
    figure(1);
    hold on;
    
    % 5HT    
    col = getCol4Stim(1, param);
    scatter(base_mango_5ht(p_mango_5ht==0), drug_mango_5ht(p_mango_5ht==0), ...
        sz, 'o', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(base_mango_5ht(p_mango_5ht==1), drug_mango_5ht(p_mango_5ht==1), ...
        sz, 'o', 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;

    
    % NaCl
    col = getCol4Stim(0, param);
    scatter(base_mango_nacl(p_mango_nacl==0), drug_mango_nacl(p_mango_nacl==0), ...
        sz, 'o', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(base_mango_nacl(p_mango_nacl==1), drug_mango_nacl(p_mango_nacl==1), ...
        sz, 'o', 'markerfacecolor', col, 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
        
    figure(2);
    hold on;
    
    % 5HT    
    col = getCol4Stim(1, param);
    scatter(rr_mango_5ht(p_mango_5ht==0), ...
        drug_mango_5ht(p_mango_5ht==0) - base_mango_5ht(p_mango_5ht==0),...
        sz, 'o', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(rr_mango_5ht(p_mango_5ht==1), ...
        drug_mango_5ht(p_mango_5ht==1) - base_mango_5ht(p_mango_5ht==1), ...
        sz, 'o', 'markerfacecolor', col, 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    
    % NaCl
    col = getCol4Stim(0, param);
    scatter(rr_mango_nacl(p_mango_nacl==0),...
        drug_mango_nacl(p_mango_nacl==0) - base_mango_nacl(p_mango_nacl==0), ...
        sz, 'o', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(rr_mango_nacl(p_mango_nacl==1),...
        drug_mango_nacl(p_mango_nacl==1) - base_mango_nacl(p_mango_nacl==1), ...
        sz, 'o', 'markerfacecolor', col, 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    
    
    % kaki =========================
    % 5HT
    kaki_5ht = find([exinfo_temp.ismango]==0 & [exinfo_temp.is5HT]==1);
    base_kaki_5ht = base(kaki_5ht);
    drug_kaki_5ht = drug(kaki_5ht);
    rr_kaki_5ht = rr(kaki_5ht);
    ss_kaki_5ht = ss(kaki_5ht);
    ss_drug_kaki_5ht = ss_drug(kaki_5ht);
    p_kaki_5ht = ones(1, length(kaki_5ht));
    p_kaki_5ht([exinfo_temp(kaki_5ht).p_mod] >= 0.05) = 0;
    disp(['Kaki, 5HT pairs, significant: ' num2str(length(find(p_kaki_5ht==1)))])
    disp(['Kaki, 5HT pairs, non-significant: ' num2str(length(find(p_kaki_5ht==0)))])
    
    base_kaki_5ht_all = [base_kaki_5ht_all base_kaki_5ht];
    drug_kaki_5ht_all = [drug_kaki_5ht_all drug_kaki_5ht];
    p_kaki_5ht_all = [p_kaki_5ht_all p_kaki_5ht];
    ss_kaki_5ht_all = [ss_kaki_5ht_all ss_kaki_5ht];
    ss_drug_kaki_5ht_all = [ss_drug_kaki_5ht_all ss_drug_kaki_5ht];
    rr_kaki_5ht_all = [rr_kaki_5ht_all, rr_kaki_5ht];
    
    % NaCl
    kaki_nacl = find([exinfo_temp.ismango]==0 & [exinfo_temp.is5HT]==0);
    base_kaki_nacl = base(kaki_nacl);
    drug_kaki_nacl = drug(kaki_nacl);
    rr_kaki_nacl = rr(kaki_nacl);
    ss_kaki_nacl = ss(kaki_nacl);
    ss_drug_kaki_nacl = ss_drug(kaki_nacl);
    p_kaki_nacl = ones(1, length(kaki_nacl));
    p_kaki_nacl([exinfo_temp(kaki_nacl).p_mod] >= 0.05) = 0;
    disp(['Kaki, NaCl pairs, significant: ' num2str(length(find(p_kaki_nacl==1)))])
    disp(['Kaki, NaCl pairs, non-significant: ' num2str(length(find(p_kaki_nacl==0)))])
    
    base_kaki_nacl_all = [base_kaki_nacl_all base_kaki_nacl];
    drug_kaki_nacl_all = [drug_kaki_nacl_all drug_kaki_nacl];
    p_kaki_nacl_all = [p_kaki_nacl_all p_kaki_nacl];
    ss_kaki_nacl_all = [ss_kaki_nacl_all ss_kaki_nacl];
    ss_drug_kaki_nacl_all = [ss_drug_kaki_nacl_all ss_drug_kaki_nacl];
    rr_kaki_nacl_all = [rr_kaki_nacl_all, rr_kaki_nacl];
    
    disp('++++++++++++++++++++++++++++')
    
    % plot -------------------------
    h(1) = figure(1);
    
    % 5HT    
    col = getCol4Stim(1, param);
    scatter(base_kaki_5ht(p_kaki_5ht==0), drug_kaki_5ht(p_kaki_5ht==0), ...
        sz, 's', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(base_kaki_5ht(p_kaki_5ht==1), drug_kaki_5ht(p_kaki_5ht==1), ...
        sz, 's', 'markerfacecolor', col, 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    
    % NaCl
    col = getCol4Stim(0, param);
    scatter(base_kaki_nacl(p_kaki_nacl==0), drug_kaki_nacl(p_kaki_nacl==0), ...
        sz, 's', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(base_kaki_nacl(p_kaki_nacl==1), drug_kaki_nacl(p_kaki_nacl==1), ...
        sz, 's', 'markerfacecolor', col, 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    
    axis([xl yl])
    
    
    h(2) = figure(2);
    % 5HT    
    col = getCol4Stim(1, param);
    scatter(rr_kaki_5ht(p_kaki_5ht==0), ...
        drug_kaki_5ht(p_kaki_5ht==0) - base_kaki_5ht(p_kaki_5ht==0),...
        sz, 's', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(rr_kaki_5ht(p_kaki_5ht==1), ...
        drug_kaki_5ht(p_kaki_5ht==1) - base_kaki_5ht(p_kaki_5ht==1), ...
        sz, 's', 'markerfacecolor', col, 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    
    % NaCl
    col = getCol4Stim(0, param);
    scatter(rr_kaki_nacl(p_kaki_nacl==0),...
        drug_kaki_nacl(p_kaki_nacl==0) - base_kaki_nacl(p_kaki_nacl==0), ...
        sz, 's', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(rr_kaki_nacl(p_kaki_nacl==1),...
        drug_kaki_nacl(p_kaki_nacl==1) - base_kaki_nacl(p_kaki_nacl==1), ...
        sz, 's', 'markerfacecolor', col, 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    
    set(gca, 'XScale','log')
    axis([xl2 yl2])
    
    % table -----------------------------
    % Mango
    [~,htp,~,htstats] = ttest(base_mango_5ht, drug_mango_5ht);
    [htpn,~,htnstats] = nansignrank(base_mango_5ht, drug_mango_5ht);
    [~,nlp,~,nlstats] = ttest(base_mango_nacl, drug_mango_nacl);
    [nlpn,~,nlnstats] = nansignrank(base_mango_nacl, drug_mango_nacl);
    [~,dp,~,dstats] = ttest2(drug_mango_5ht - base_mango_5ht, drug_mango_nacl - base_mango_nacl);
    [dpn,~,dnstats] = nanranksum(drug_mango_5ht - base_mango_5ht, drug_mango_nacl - base_mango_nacl);
    [rho01, rhop01] = nancorr(base_mango_5ht', drug_mango_5ht', 'Spearman');
    [rho11, rhop11] = nancorr([base_mango_5ht - drug_mango_5ht]', rr_mango_5ht', 'Spearman');
    [rho02, rhop02] = nancorr(base_mango_nacl', drug_mango_nacl', 'Spearman');
    [rho12, rhop12] = nancorr([base_mango_nacl - drug_mango_nacl]', rr_mango_nacl', 'Spearman');

    T(:,i+4) = [nanlength(base_mango_5ht) + nanlength(base_mango_nacl), ...
        nanlength(base_mango_5ht), nanlength(base_mango_nacl), ...
        nonzeromean([ss_mango_5ht ss_drug_mango_5ht ss_mango_nacl ss_drug_mango_nacl]),...
        nonzeromedian([ss_mango_5ht ss_drug_mango_5ht ss_mango_nacl ss_drug_mango_nacl]),...
        nanmax(ss_mango_5ht), nanmax(ss_drug_mango_5ht),...
        nanmax(ss_mango_nacl), nanmax(ss_drug_mango_nacl),...
        nanmin(ss_mango_5ht), nanmin(ss_drug_mango_5ht),...
        nanmin(ss_mango_nacl), nanmin(ss_drug_mango_nacl),...
        nanmedian([base_mango_5ht base_mango_nacl]),...
        nanmedian([drug_mango_5ht drug_mango_nacl]),...
        rho01, rhop01, rho11, rhop11,...
        htstats.tstat, htp, htnstats.signedrank, htpn, ...     
        length(find(p_mango_5ht==1 & ~isnan(base_mango_5ht))),length(find(p_mango_5ht==0 & ~isnan(base_mango_5ht))),...
        nanmedian(base_mango_5ht), nanmedian(drug_mango_5ht),...
        rho02, rhop02, rho12, rhop12,...
        nlstats.tstat, nlp, nlnstats.signedrank, nlpn,...
        length(find(p_mango_nacl==1 & ~isnan(base_mango_nacl))),length(find(p_mango_nacl==0 & ~isnan(base_mango_nacl))),...
        nanmedian(base_mango_nacl), nanmedian(drug_mango_nacl),...
        dstats.tstat, dp, dnstats.ranksum, dpn];
    
    % Kaki
    [~,htp,~,htstats] = ttest(base_kaki_5ht, drug_kaki_5ht);
    [htpn,~,htnstats] = nansignrank(base_kaki_5ht, drug_kaki_5ht);
    [~,nlp,~,nlstats] = ttest(base_kaki_nacl, drug_kaki_nacl);
    [nlpn,~,nlnstats] = nansignrank(base_kaki_nacl, drug_kaki_nacl);
    [~,dp,~,dstats] = ttest2(drug_kaki_5ht - base_kaki_5ht, drug_kaki_nacl - base_kaki_nacl);
    [dpn,~,dnstats] = nanranksum(drug_kaki_5ht - base_kaki_5ht, drug_kaki_nacl - base_kaki_nacl);
    [rho01, rhop01] = nancorr(base_kaki_5ht', drug_kaki_5ht', 'Spearman');
    [rho11, rhop11] = nancorr([base_kaki_5ht - drug_kaki_5ht]', rr_kaki_5ht', 'Spearman');
    [rho02, rhop02] = nancorr(base_kaki_nacl', drug_kaki_nacl', 'Spearman');
    [rho12, rhop12] = nancorr([base_kaki_nacl - drug_kaki_nacl]', rr_kaki_nacl', 'Spearman');
    
    T(:,i+8) = [nanlength(base_kaki_5ht) + nanlength(base_kaki_nacl), ...
        nanlength(base_kaki_5ht), nanlength(base_kaki_nacl), ...
        nonzeromean([ss_kaki_5ht ss_drug_kaki_5ht ss_kaki_nacl ss_drug_kaki_nacl]),...
        nonzeromedian([ss_kaki_5ht ss_drug_kaki_5ht ss_kaki_nacl ss_drug_kaki_nacl]),...
        nanmax(ss_kaki_5ht), nanmax(ss_drug_kaki_5ht),...
        nanmax(ss_kaki_nacl), nanmax(ss_drug_kaki_nacl),...
        nanmin(ss_kaki_5ht), nanmin(ss_drug_kaki_5ht),...
        nanmin(ss_kaki_nacl), nanmin(ss_drug_kaki_nacl),...
        nanmedian([base_kaki_5ht base_kaki_nacl]),...
        nanmedian([drug_kaki_5ht drug_kaki_nacl]),...
        rho01, rhop01, rho11, rhop11,...
        htstats.tstat, htp, htnstats.signedrank, htpn, ...     
        length(find(p_kaki_5ht==1 & ~isnan(base_kaki_5ht))),length(find(p_kaki_5ht==0 & ~isnan(base_kaki_5ht))),...
        nanmedian(base_kaki_5ht), nanmedian(drug_kaki_5ht),...
        rho02, rhop02, rho12, rhop12,...
        nlstats.tstat, nlp, nlnstats.signedrank, nlpn,...
        length(find(p_kaki_nacl==1 & ~isnan(base_kaki_nacl))),length(find(p_kaki_nacl==0 & ~isnan(base_kaki_nacl))),...
        nanmedian(base_kaki_nacl), nanmedian(drug_kaki_nacl),...
        dstats.tstat, dp, dnstats.ranksum, dpn];
    
    % Both
    [~,htp,~,htstats] = ttest([base_mango_5ht base_kaki_5ht], [drug_mango_5ht drug_kaki_5ht]);
    [htpn,~,htnstats] = nansignrank([base_mango_5ht base_kaki_5ht], [drug_mango_5ht drug_kaki_5ht]);
    [~,nlp,~,nlstats] = ttest([base_mango_nacl base_kaki_nacl], [drug_mango_nacl drug_kaki_nacl]);
    [nlpn,~,nlnstats] = nansignrank([base_mango_nacl base_kaki_nacl], [drug_mango_nacl drug_kaki_nacl]);
    [~,dp,~,dstats] = ttest2([drug_mango_5ht - base_mango_5ht drug_kaki_5ht - base_kaki_5ht],...
     [drug_mango_nacl - base_mango_nacl drug_kaki_nacl - base_kaki_nacl]);
     [dpn,~,dnstats] = nanranksum([drug_mango_5ht - base_mango_5ht drug_kaki_5ht - base_kaki_5ht],...
     [drug_mango_nacl - base_mango_nacl drug_kaki_nacl - base_kaki_nacl]);
    [rho01, rhop01] = nancorr([base_mango_5ht base_kaki_5ht]', [drug_mango_5ht drug_kaki_5ht]', 'Spearman');
    [rho11, rhop11] = nancorr([base_mango_5ht - drug_mango_5ht base_kaki_5ht - drug_kaki_5ht]', [rr_mango_5ht rr_kaki_5ht]', 'Spearman');
    [rho02, rhop02] = nancorr([base_mango_nacl base_kaki_nacl]', [drug_mango_nacl drug_kaki_nacl]', 'Spearman');
    [rho12, rhop12] = nancorr([base_mango_nacl - drug_mango_nacl base_kaki_nacl - drug_kaki_nacl]', [rr_mango_nacl rr_kaki_nacl]', 'Spearman');
    
     T(:,i) = [nanlength(base_mango_5ht) + nanlength(base_mango_nacl) + nanlength(base_kaki_5ht) + nanlength(base_kaki_nacl), ...
        nanlength(base_mango_5ht) + nanlength(base_kaki_5ht), nanlength(base_mango_nacl) + nanlength(base_kaki_nacl), ...
        nonzeromean([ss_mango_5ht ss_drug_mango_5ht ss_mango_nacl ss_drug_mango_nacl ...
        ss_kaki_5ht ss_drug_kaki_5ht ss_kaki_nacl ss_drug_kaki_nacl]),...
        nonzeromedian([ss_mango_5ht ss_drug_mango_5ht ss_mango_nacl ss_drug_mango_nacl ...
        ss_kaki_5ht ss_drug_kaki_5ht ss_kaki_nacl ss_drug_kaki_nacl]),...
        nanmax([ss_mango_5ht ss_kaki_5ht]), nanmax([ss_drug_mango_5ht ss_drug_kaki_5ht]),...
        nanmax([ss_mango_nacl ss_kaki_nacl]), nanmax([ss_drug_mango_nacl ss_drug_kaki_nacl]),...
        nanmin([ss_mango_5ht ss_kaki_5ht]), nanmin([ss_drug_mango_5ht ss_drug_kaki_5ht]),...
        nanmin([ss_mango_nacl ss_kaki_nacl]), nanmin([ss_drug_mango_nacl ss_drug_kaki_nacl]),...
        nanmedian([base_mango_5ht base_mango_nacl base_kaki_5ht base_kaki_nacl]),...
        nanmedian([drug_mango_5ht drug_mango_nacl drug_kaki_5ht drug_kaki_nacl]),...
        rho01, rhop01, rho11, rhop11,...
        htstats.tstat, htp, htnstats.signedrank, htpn, ...     
        length(find(p_mango_5ht==1 & ~isnan(base_mango_5ht))) + length(find(p_kaki_5ht==1 & ~isnan(base_kaki_5ht))),...
        length(find(p_mango_5ht==0 & ~isnan(drug_mango_5ht))) + length(find(p_kaki_5ht==0 & ~isnan(drug_kaki_5ht))),...
        nanmedian([base_mango_5ht base_kaki_5ht]), nanmedian([drug_mango_5ht drug_kaki_5ht]),...
        rho02, rhop02, rho12, rhop12,...
        nlstats.tstat, nlp, nlnstats.signedrank, nlpn,...
        length(find(p_mango_nacl==1 & ~isnan(base_mango_nacl))) + length(find(p_kaki_nacl==1 & ~isnan(base_kaki_nacl))),...
        length(find(p_mango_nacl==0 & ~isnan(base_mango_nacl))) + length(find(p_kaki_nacl==0 & ~isnan(base_kaki_nacl))),...
        nanmedian([base_mango_nacl base_kaki_nacl]), nanmedian([drug_mango_nacl drug_kaki_nacl]),...
        dstats.tstat, dp, dnstats.ranksum, dpn];
    
end

% store stimulus collapsed data -----------------------
% Mango
[~,htp,~,htstats] = ttest(base_mango_5ht_all, drug_mango_5ht_all);
[htpn,~,htnstats] = nansignrank(base_mango_5ht_all, drug_mango_5ht_all);
[~,nlp,~,nlstats] = ttest(base_mango_nacl_all, drug_mango_nacl_all);
[nlpn,~,nlnstats] = nansignrank(base_mango_nacl_all, drug_mango_nacl_all);
[~,dp,~,dstats] = ttest2(drug_mango_5ht_all - base_mango_5ht_all, drug_mango_nacl_all - base_mango_nacl_all);
[dpn,~,dnstats] = nanranksum(drug_mango_5ht_all - base_mango_5ht_all, drug_mango_nacl_all - base_mango_nacl_all);
[rho01, rhop01] = nancorr(base_mango_5ht_all', drug_mango_5ht_all', 'Spearman');
[rho11, rhop11] = nancorr([base_mango_5ht_all - drug_mango_5ht_all]', rr_mango_5ht_all', 'Spearman');
[rho02, rhop02] = nancorr(base_mango_nacl_all', drug_mango_nacl_all', 'Spearman');
[rho12, rhop12] = nancorr([base_mango_nacl_all - drug_mango_nacl_all]', rr_mango_nacl_all', 'Spearman');

T(:,14) = [nanlength([base_mango_5ht_all base_mango_nacl_all]), ...
    nanlength(base_mango_5ht_all), nanlength(base_mango_nacl_all),...
    nonzeromean([ss_mango_5ht_all ss_drug_mango_5ht_all ss_mango_nacl_all ss_drug_mango_nacl_all]),...    
    nonzeromedian([ss_mango_5ht_all ss_drug_mango_5ht_all ss_mango_nacl_all ss_drug_mango_nacl_all]),...
    nanmax(ss_mango_5ht_all),nanmax(ss_drug_mango_5ht_all),...
    nanmax(ss_mango_nacl_all),nanmax(ss_drug_mango_nacl_all),...
    nanmin(ss_mango_5ht_all),nanmin(ss_drug_mango_5ht_all),...
    nanmin(ss_mango_nacl_all),nanmin(ss_drug_mango_nacl_all),...
    nanmedian([base_mango_5ht_all base_mango_nacl_all]),...    
    nanmedian([drug_mango_5ht_all drug_mango_nacl_all]),...
    rho01, rhop01, rho11, rhop11,...
    htstats.tstat, htp, htnstats.signedrank, htpn, ...     
    length(find(p_mango_5ht_all==1 & ~isnan(base_mango_5ht_all))), ...
    length(find(p_mango_5ht_all==0 & ~isnan(base_mango_5ht_all))), ...
    nanmedian(base_mango_5ht_all), nanmedian(drug_mango_5ht_all),...
    rho02, rhop02, rho12, rhop12,...
    nlstats.tstat, nlp, nlnstats.signedrank, nlpn,...
    length(find(p_mango_nacl_all==1 & ~isnan(base_mango_nacl_all))), ...
    length(find(p_mango_nacl_all==0 & ~isnan(base_mango_nacl_all))),...
    nanmedian(base_mango_nacl_all), nanmedian(drug_mango_nacl_all),...
    dstats.tstat, dp, dnstats.ranksum, dpn];

% Kaki
[~,htp,~,htstats] = ttest(base_kaki_5ht_all, drug_kaki_5ht_all);
[htpn,~,htnstats] = nansignrank(base_kaki_5ht_all, drug_kaki_5ht_all);
[~,nlp,~,nlstats] = ttest(base_kaki_nacl_all, drug_kaki_nacl_all);
[nlpn,~,nlnstats] = nansignrank(base_kaki_nacl_all, drug_kaki_nacl_all);
[~,dp,~,dstats] = ttest2(drug_kaki_5ht_all - base_kaki_5ht_all, drug_kaki_nacl_all - base_kaki_nacl_all);
[dpn,~,dnstats] = nanranksum(drug_kaki_5ht_all - base_kaki_5ht_all, drug_kaki_nacl_all - base_kaki_nacl_all);
[rho01, rhop01] = nancorr(base_kaki_5ht_all', drug_kaki_5ht_all', 'Spearman');
[rho11, rhop11] = nancorr([base_kaki_5ht_all - drug_kaki_5ht_all]', rr_kaki_5ht_all', 'Spearman');
[rho02, rhop02] = nancorr(base_kaki_nacl_all', drug_kaki_nacl_all', 'Spearman');
[rho12, rhop12] = nancorr([base_kaki_nacl_all - drug_kaki_nacl_all]', rr_kaki_nacl_all', 'Spearman');

T(:,15) = [nanlength([base_kaki_5ht_all base_kaki_nacl_all]), ...
    nanlength(base_kaki_5ht_all), nanlength(base_kaki_nacl_all),...
    nonzeromean([ss_kaki_5ht_all ss_drug_kaki_5ht_all ss_kaki_nacl_all ss_drug_kaki_nacl_all]),...    
    nonzeromedian([ss_kaki_5ht_all ss_drug_kaki_5ht_all ss_kaki_nacl_all ss_drug_kaki_nacl_all]),...
    nanmax(ss_kaki_5ht_all),nanmax(ss_drug_kaki_5ht_all),...
    nanmax(ss_kaki_nacl_all),nanmax(ss_drug_kaki_nacl_all),...
    nanmin(ss_kaki_5ht_all),nanmin(ss_drug_kaki_5ht_all),...
    nanmin(ss_kaki_nacl_all),nanmin(ss_drug_kaki_nacl_all),...
    nanmedian([base_kaki_5ht_all base_kaki_nacl_all]),...    
    nanmedian([drug_kaki_5ht_all drug_kaki_nacl_all]),...
    rho01, rhop01, rho11, rhop11,...
    htstats.tstat, htp, htnstats.signedrank, htpn, ...     
    length(find(p_kaki_5ht_all==1 & ~isnan(base_kaki_5ht_all))), ...
    length(find(p_kaki_5ht_all==0 & ~isnan(base_kaki_5ht_all))), ...
    nanmedian(base_kaki_5ht_all), nanmedian(drug_kaki_5ht_all),...
    rho02, rhop02, rho12, rhop12,...
    nlstats.tstat, nlp, nlnstats.signedrank, nlpn,...
    length(find(p_kaki_nacl_all==1 & ~isnan(base_kaki_nacl_all))), ...
    length(find(p_kaki_nacl_all==0 & ~isnan(base_kaki_nacl_all))),...
    nanmedian(base_kaki_nacl_all), nanmedian(drug_kaki_nacl_all),...
    dstats.tstat, dp, dnstats.ranksum, dpn];

% Both
[~,htp,~,htstats] = ttest([base_mango_5ht_all base_kaki_5ht_all], [drug_mango_5ht_all drug_kaki_5ht_all]);
[htpn,~,htnstats] = nansignrank([base_mango_5ht_all base_kaki_5ht_all], [drug_mango_5ht_all drug_kaki_5ht_all]);
[~,nlp,~,nlstats] = ttest([base_mango_nacl_all base_kaki_nacl_all], [drug_mango_nacl_all drug_kaki_nacl_all]);
[nlpn,~,nlnstats] = nansignrank([base_mango_nacl_all base_kaki_nacl_all], [drug_mango_nacl_all drug_kaki_nacl_all]);
[~,dp,~,dstats] = ttest2([drug_mango_5ht_all - base_mango_5ht_all drug_kaki_5ht_all - base_kaki_5ht_all], ...
    [drug_mango_nacl_all - base_mango_nacl_all drug_kaki_nacl_all - base_kaki_nacl_all]);
[dpn,~,dnstats] = nanranksum([drug_mango_5ht_all - base_mango_5ht_all drug_kaki_5ht_all - base_kaki_5ht_all], ...
     [drug_mango_nacl_all - base_mango_nacl_all drug_kaki_nacl_all - base_kaki_nacl_all]);
[rho01, rhop01] = nancorr([base_mango_5ht_all base_kaki_5ht_all]', [drug_mango_5ht_all drug_kaki_5ht_all]', 'Spearman');
[rho11, rhop11] = nancorr([base_mango_5ht_all - drug_mango_5ht_all base_kaki_5ht_all - drug_kaki_5ht_all]', [rr_mango_5ht_all rr_kaki_5ht_all]', 'Spearman');
[rho02, rhop02] = nancorr([base_mango_nacl_all base_kaki_nacl_all]', [drug_mango_nacl_all drug_kaki_nacl_all]', 'Spearman');
[rho12, rhop12] = nancorr([base_mango_nacl_all - drug_mango_nacl_all base_kaki_nacl_all - drug_kaki_nacl_all]', [rr_mango_nacl_all rr_kaki_nacl_all]', 'Spearman');

 T(:,13) = [nanlength([base_mango_5ht_all base_mango_nacl_all]) + nanlength([base_kaki_5ht_all base_kaki_nacl_all]), ...
    nanlength(base_mango_5ht_all)+nanlength(base_kaki_5ht_all), nanlength(base_mango_nacl_all) + nanlength(base_kaki_nacl_all),...
    nonzeromean([ss_mango_5ht_all ss_drug_mango_5ht_all ss_mango_nacl_all ss_drug_mango_nacl_all, ...
    ss_kaki_5ht_all ss_drug_kaki_5ht_all ss_kaki_nacl_all ss_drug_kaki_nacl_all]),...    
    nonzeromedian([ss_mango_5ht_all ss_drug_mango_5ht_all ss_mango_nacl_all ss_drug_mango_nacl_all, ...
    ss_kaki_5ht_all ss_drug_kaki_5ht_all ss_kaki_nacl_all ss_drug_kaki_nacl_all]),...    
    nanmax([ss_mango_5ht_all ss_kaki_5ht_all]),nanmax([ss_drug_mango_5ht_all ss_drug_kaki_5ht_all]),...
    nanmax([ss_mango_nacl_all ss_kaki_nacl_all]),nanmax([ss_drug_mango_nacl_all ss_drug_kaki_nacl_all]),...
    nanmin([ss_mango_5ht_all ss_kaki_5ht_all]),nanmin([ss_drug_mango_5ht_all ss_drug_kaki_5ht_all]),...
    nanmin([ss_mango_nacl_all ss_kaki_nacl_all]),nanmin([ss_drug_mango_nacl_all ss_drug_kaki_nacl_all]),...
    nanmedian([base_mango_5ht_all base_mango_nacl_all base_kaki_5ht_all base_kaki_nacl_all]),...    
    nanmedian([drug_mango_5ht_all drug_mango_nacl_all drug_kaki_5ht_all drug_kaki_nacl_all]),...
    rho01, rhop01, rho11, rhop11,...
    htstats.tstat, htp, htnstats.signedrank, htpn, ...     
    length(find(p_mango_5ht_all==1 & ~isnan(base_mango_5ht_all))) + length(find(p_kaki_5ht_all==1 & ~isnan(base_kaki_5ht_all))), ...
    length(find(p_mango_5ht_all==0 & ~isnan(base_mango_5ht_all))) + length(find(p_kaki_5ht_all==0 & ~isnan(base_kaki_5ht_all))), ...
    nanmedian([base_mango_5ht_all base_kaki_5ht_all]), nanmedian([drug_mango_5ht_all drug_kaki_5ht_all]),...
    rho02, rhop02, rho12, rhop12,...
    nlstats.tstat, nlp, nlnstats.signedrank, nlpn,...
    length(find(p_mango_nacl_all==1 & ~isnan(base_mango_nacl_all))) + length(find(p_kaki_nacl_all==1 & ~isnan(base_kaki_nacl_all))), ...
    length(find(p_mango_nacl_all==0 & ~isnan(base_mango_nacl_all))) + length(find(p_kaki_nacl_all==0 & ~isnan(base_kaki_nacl_all))),...
    nanmedian([base_mango_nacl_all base_kaki_nacl_all]), nanmedian([drug_mango_nacl_all drug_kaki_nacl_all]),...
    dstats.tstat, dp, dnstats.ranksum, dpn];

% sample size summary
% per experiment
disp([name '_' num2str(nanminss) '+++++++'])
disp('--------------------------------------------')
disp('sample size per experiment:')
disp(['median: ' num2str(nonzeromedian(ssperexp))])
disp(['max: ' num2str(nanmax(ssperexp))])
disp(['min: ' num2str(nanmin(ssperexp))])

% per stimulus
disp('--------------------------------------------')
disp('sample size per stimulus:')
disp(['median: ' num2str(nonzeromedian(ssperstm))])
disp(['max: ' num2str(nanmax(ssperstm))])
disp(['min: ' num2str(nanmin(ssperstm))])

% % correlation to relative ratio
% ser = [base_mango_5ht_all - drug_mango_5ht_all, base_kaki_5ht_all - drug_kaki_5ht_all];
% ser(isnan(ser)) = [];
% list_inc(isnan(list_inc)) = [];
% rr_ser = [exinfo(intersect(find([exinfo.is5HT]==1), list_inc)).nonparam_ratio];
% % rr_ser = [rr_mango_5ht_all rr_kaki_5ht_all];
% [rho, pp] = corr(ser', rr_ser', 'type', 'Spearman');
% disp(['correlation with ' name ' (5HT) with the relative ratio:'])
% disp(['rho = ' num2str(rho) ', p = ' num2str(pp)])

if saveoption==1
    
    % save figure
    savefig(h(1), ['Z:\Corinna\SerotoninPaper\Figures\Figure09_NoiseCorr_and_FF\raw_figs\' name '_' num2str(nanminss) '.fig'])
    savefig(h(2), ['Z:\Corinna\SerotoninPaper\Figures\Figure09_NoiseCorr_and_FF\raw_figs\' name '_vsRelativeGain_' num2str(nanminss) '.fig'])
    disp(['The figure was saved in '  ['Z:\Corinna\SerotoninPaper\Figures\Figure09_NoiseCorr_and_FF\raw_figs\' name '_' num2str(nanminss)]])
    
    % save table
    VariableNames = {'or_both','sf_both','co_both', 'sz_both',...
        'or_mango','sf_mango', 'co_mango', 'sz_mango',...
        'or_kaki','sf_kaki', 'co_kaki', 'sz_kaki','collapsed_all','mango_all','kaki_all'} ;
    RowNames = {'total', 'n_5HT', 'n_NaCl','mean_ss','median_ss',...
        'nanmax_SS_Base_5ht','nanmax_SS_5HT','nanmax_SS_Base_nacl','nanmax_SS_NaCl',...
        'nanmin_SS_Base_5ht','nanmin_SS_5HT','nanmin_SS_Base_nacl','nanmin_SS_NaCl',...
        'base; median','drug; median',...
        '5HT: rankcorr, r = ', '5HT: rankcorr, p = ','5HTvsRR: rankcorr, r = ', '5HTvsRR: rankcorr, p = ',...
        '5HT: ttest, t = ', '5HT: ttest, p = ', '5HT: signrank, rank = ','5HT: signrank, p = ', ...
        '5ht; significant', '5ht; non-significant', '5ht_base; median','5ht; median', ...
        'nacl: rankcorr, r = ', 'nacl: rankcorr, p = ','naclvsRR: rankcorr, r = ', 'naclvsRR: rankcorr, p = ', ...
        'nacl: ttest, t = ', 'nacl: ttest, p = ', 'nacl: signrank, rank = ','nacl: signrank, p = ', ...
        'nacl; significant', 'nacl; non-significant',   'nacl_base; median','nacl; median',...
        'diff: ttest2, t = ', 'diff: ttest, p = ', 'diff: ranksum, rank = ','diff: ranksum, p = '};
    T = array2table(T, 'RowNames', RowNames,'VariableNames',VariableNames);
    writetable(T, ['Z:\Corinna\SerotoninPaper\Figures\Figure09_NoiseCorr_and_FF\raw_figs\' name '_' num2str(nanminss) '.csv'],'WriteRowNames',true)
end
