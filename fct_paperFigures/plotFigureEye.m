function h = plotFigureEye(exinfo, name, saveoption, varargin)
%% plot scatter for eye data for revision
% INPUT: exinfo ... given by Corinna's 'runExinfoAnalysis'
%             name ... 'fixationspan' or 'microsac_amplitude' or
%             'microsac_counts'
%
% written by Katsuhisa (13.07.17)
% ++++++++++++++++++++++++++++++++++++

close all

switch nargin
    case 0
        error('provide exinfo!')
    case 1
        name = 'fixation_precision';
        saveoption = 0;
    case 2
        saveoption = 0;
end

if strcmp(name, 'fixation_precision')
    name = 'fixationspan';
end

% basic figure parameters
sz = 10;
a = 0.5;

switch name
    case 'fixationspan'
        xl = [0.2 1];
        yl = [0.2 1];
    case 'microsac_amplitude'
        xl = [0.03 0.12];
        yl = [0.03 0.12];
    case 'microsac_counts'
        xl = [0.5 25];
        yl = [0.5 25];
end
 

% load list of index
% load('Z:\Corinna\SharedCode\Katsu\incl_i_one_stimulus_cond.mat')
load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stim_cond_2007.mat')
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
T = nan(25,15);
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
h = setPaperPlotsProps();
plot(xl, yl, '-', 'color', 0.6*ones(1,3),'linewidth',0.5)
hold on;
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
    
    exinfo_temp = exinfo([exinfo.stmtype]==stmtype);
    switch name
        case 'fixationspan'
            base = arrayfun(@(x) x.(name)(1), exinfo_temp);
            drug = arrayfun(@(x) x.([name '_drug'])(1), exinfo_temp);
        case 'microsac_amplitude'
            base = nan(1, length(exinfo_temp));
            drug = nan(1, length(exinfo_temp));
            for k = 1:length(exinfo_temp)
                base(k) = median(exinfo_temp(k).(name));
                drug(k) = median(exinfo_temp(k).([name '_drug']));
            end
        case 'microsac_counts'
            base = nan(1, length(exinfo_temp));
            drug = nan(1, length(exinfo_temp));
            for k = 1:length(exinfo_temp)
%                 base(k) = log(mean(exinfo_temp(k).(name)));
%                 drug(k) = log(mean(exinfo_temp(k).([name '_drug'])));
                base(k) = mean(exinfo_temp(k).(name));
                drug(k) = mean(exinfo_temp(k).([name '_drug']));
            end
            set(gca,'xscale','log')
            set(gca,'yscale','log')
    end
    
    % mango =======================
    % 5HT
    mango_5ht = find([exinfo_temp.ismango]==1 & [exinfo_temp.is5HT]==1);
    base_mango_5ht = base(mango_5ht);
    drug_mango_5ht = drug(mango_5ht);
    p_mango_5ht = ones(1, length(mango_5ht));
    p_mango_5ht([exinfo_temp(mango_5ht).p_mod] >= 0.05) = 0;
    disp(['Mango, 5HT pairs, significant: ' num2str(length(find(p_mango_5ht==1)))])
    disp(['Mango, 5HT pairs, non-significant: ' num2str(length(find(p_mango_5ht==0)))])
    
    base_mango_5ht_all = [base_mango_5ht_all base_mango_5ht];
    drug_mango_5ht_all = [drug_mango_5ht_all drug_mango_5ht];
    p_mango_5ht_all = [p_mango_5ht_all p_mango_5ht];
    
    % NaCl
    mango_nacl = find([exinfo_temp.ismango]==1 & [exinfo_temp.is5HT]==0);
    base_mango_nacl = base(mango_nacl);
    drug_mango_nacl = drug(mango_nacl);
    p_mango_nacl = ones(1, length(mango_nacl));
    p_mango_nacl([exinfo_temp(mango_nacl).p_mod] >= 0.05) = 0;
    disp(['Mango, NaCl pairs, significant: ' num2str(length(find(p_mango_nacl==1)))])
    disp(['Mango, NaCl pairs, non-significant: ' num2str(length(find(p_mango_nacl==0)))])
    
    base_mango_nacl_all = [base_mango_nacl_all base_mango_nacl];
    drug_mango_nacl_all = [drug_mango_nacl_all drug_mango_nacl];
    p_mango_nacl_all = [p_mango_nacl_all p_mango_nacl];
    
    % plot -------------------------
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
        sz, 'o', 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
        
    
    % kaki =========================
    % 5HT
    kaki_5ht = find([exinfo_temp.ismango]==0 & [exinfo_temp.is5HT]==1);
    base_kaki_5ht = base(kaki_5ht);
    drug_kaki_5ht = drug(kaki_5ht);
    p_kaki_5ht = ones(1, length(kaki_5ht));
    p_kaki_5ht([exinfo_temp(kaki_5ht).p_mod] >= 0.05) = 0;
    disp(['Kaki, 5HT pairs, significant: ' num2str(length(find(p_kaki_5ht==1)))])
    disp(['Kaki, 5HT pairs, non-significant: ' num2str(length(find(p_kaki_5ht==0)))])
    
    base_kaki_5ht_all = [base_kaki_5ht_all base_kaki_5ht];
    drug_kaki_5ht_all = [drug_kaki_5ht_all drug_kaki_5ht];
    p_kaki_5ht_all = [p_kaki_5ht_all p_kaki_5ht];
    
    % NaCl
    kaki_nacl = find([exinfo_temp.ismango]==0 & [exinfo_temp.is5HT]==0);
    base_kaki_nacl = base(kaki_nacl);
    drug_kaki_nacl = drug(kaki_nacl);
    p_kaki_nacl = ones(1, length(kaki_nacl));
    p_kaki_nacl([exinfo_temp(kaki_nacl).p_mod] >= 0.05) = 0;
    disp(['Kaki, NaCl pairs, significant: ' num2str(length(find(p_kaki_nacl==1)))])
    disp(['Kaki, NaCl pairs, non-significant: ' num2str(length(find(p_kaki_nacl==0)))])
    
    base_kaki_nacl_all = [base_kaki_nacl_all base_kaki_nacl];
    drug_kaki_nacl_all = [drug_kaki_nacl_all drug_kaki_nacl];
    p_kaki_nacl_all = [p_kaki_nacl_all p_kaki_nacl];
    
    disp('++++++++++++++++++++++++++++')
    
    % plot -------------------------
    % 5HT    
    col = getCol4Stim(1, param);
    scatter(base_kaki_5ht(p_kaki_5ht==0), drug_kaki_5ht(p_kaki_5ht==0), ...
        sz, 's', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(base_kaki_5ht(p_kaki_5ht==1), drug_kaki_5ht(p_kaki_5ht==1), ...
        sz, 's', 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    
    % NaCl
    col = getCol4Stim(0, param);
    scatter(base_kaki_nacl(p_kaki_nacl==0), drug_kaki_nacl(p_kaki_nacl==0), ...
        sz, 's', 'markerfacecolor', 'w', 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    scatter(base_kaki_nacl(p_kaki_nacl==1), drug_kaki_nacl(p_kaki_nacl==1), ...
        sz, 's', 'markerfacecolor', col, 'markeredgecolor', col, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    
    % table -----------------------------
    % Mango
    [~,htp,~,htstats] = ttest(base_mango_5ht, drug_mango_5ht);
    [htpn,~,htnstats] = nansignrank(base_mango_5ht, drug_mango_5ht);
    [~,nlp,~,nlstats] = ttest(base_mango_nacl, drug_mango_nacl);
    [nlpn,~,nlnstats] = nansignrank(base_mango_nacl, drug_mango_nacl);
    [~,dp,~,dstats] = ttest2(drug_mango_5ht - base_mango_5ht, drug_mango_nacl - base_mango_nacl);
    [dpn,~,dnstats] = nanranksum(drug_mango_5ht - base_mango_5ht, drug_mango_nacl - base_mango_nacl);
    
    T(:,i+4) = [length(mango_5ht)+length(mango_nacl), length(p_mango_5ht), length(p_mango_nacl), median([base_mango_5ht base_mango_nacl]),...
        median([drug_mango_5ht drug_mango_nacl]), htstats.tstat, htp, htnstats.signedrank, htpn, ...     
        length(find(p_mango_5ht==1)),length(find(p_mango_5ht==0)),...
        median(base_mango_5ht), median(drug_mango_5ht),...
        nlstats.tstat, nlp, nlnstats.signedrank, nlpn,  length(find(p_mango_nacl==1)),length(find(p_mango_nacl==0)),...
        median(base_mango_nacl), median(drug_mango_nacl),...
        dstats.tstat, dp, dnstats.ranksum, dpn];
    
    % Kaki
    [~,htp,~,htstats] = ttest(base_kaki_5ht, drug_kaki_5ht);
    [htpn,~,htnstats] = nansignrank(base_kaki_5ht, drug_kaki_5ht);
    [~,nlp,~,nlstats] = ttest(base_kaki_nacl, drug_kaki_nacl);
    [nlpn,~,nlnstats] = nansignrank(base_kaki_nacl, drug_kaki_nacl);
    [~,dp,~,dstats] = ttest2(drug_kaki_5ht - base_kaki_5ht, drug_kaki_nacl - base_kaki_nacl);
    [dpn,~,dnstats] = nanranksum(drug_kaki_5ht - base_kaki_5ht, drug_kaki_nacl - base_kaki_nacl);
    
    T(:,i+8) = [length(kaki_5ht)+length(kaki_nacl), length(p_kaki_5ht), length(p_kaki_nacl),median([base_kaki_5ht base_kaki_nacl]),...
        median([drug_kaki_5ht drug_kaki_nacl]), htstats.tstat, htp, htnstats.signedrank, htpn, ...     
        length(find(p_kaki_5ht==1)),length(find(p_kaki_5ht==0)),...
        median(base_kaki_5ht), median(drug_kaki_5ht),...
        nlstats.tstat, nlp, nlnstats.signedrank, nlpn,  length(find(p_kaki_nacl==1)),length(find(p_kaki_nacl==0)),...
        median(base_kaki_nacl), median(drug_kaki_nacl),...
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
    
    T(:,i) = [length(mango_5ht)+length(mango_nacl) + length(kaki_5ht)+length(kaki_nacl), ...
        length(p_mango_5ht)+length(p_kaki_5ht), length(p_mango_nacl)+length(p_kaki_nacl), ...
        median([base_mango_5ht base_mango_nacl base_kaki_5ht base_kaki_nacl]),...
        median([base_mango_5ht base_mango_nacl drug_kaki_5ht drug_kaki_nacl]),...
        htstats.tstat, htp, htnstats.signedrank, htpn, length(find([p_mango_5ht p_kaki_5ht]==1)),length(find([p_mango_5ht p_kaki_5ht]==0)),...
        median([base_mango_5ht base_kaki_5ht]), median([drug_mango_5ht drug_kaki_5ht]),...
        nlstats.tstat, nlp, nlnstats.signedrank, nlpn,  length(find([p_mango_nacl p_kaki_nacl]==1)),length(find([p_mango_nacl p_kaki_nacl]==0)),...
        median([base_mango_nacl base_kaki_nacl]), median([drug_mango_nacl drug_kaki_nacl]),...
        dstats.tstat, dp, dnstats.ranksum, dpn];    
    
end

axis([xl yl])

% store stimulus collapsed data -----------------------
% Mango
[~,htp,~,htstats] = ttest(base_mango_5ht_all, drug_mango_5ht_all);
[htpn,~,htnstats] = nansignrank(base_mango_5ht_all, drug_mango_5ht_all);
[~,nlp,~,nlstats] = ttest(base_mango_nacl_all, drug_mango_nacl_all);
[nlpn,~,nlnstats] = nansignrank(base_mango_nacl_all, drug_mango_nacl_all);
[~,dp,~,dstats] = ttest2(drug_mango_5ht_all - base_mango_5ht_all, drug_mango_nacl_all - base_mango_nacl_all);
[dpn,~,dnstats] = nanranksum(drug_mango_5ht_all - base_mango_5ht_all, drug_mango_nacl_all - base_mango_nacl_all);

T(:,14) = [length(find([exinfo.ismango]==1)), length(p_mango_5ht_all), length(p_mango_nacl_all),...
    median([base_mango_5ht_all base_mango_nacl_all]),...    
    median([drug_mango_5ht_all drug_mango_nacl_all]), htstats.tstat, htp, htnstats.signedrank, htpn, ...     
    length(find(p_mango_5ht_all==1)),length(find(p_mango_5ht_all==0)),...
    median(base_mango_5ht_all), median(drug_mango_5ht_all),...
    nlstats.tstat, nlp, nlnstats.signedrank, nlpn,  length(find(p_mango_nacl_all==1)),length(find(p_mango_nacl_all==0)),...
    median(base_mango_nacl_all), median(drug_mango_nacl_all),...
    dstats.tstat, dp, dnstats.ranksum, dpn];

% Kaki
[~,htp,~,htstats] = ttest(base_kaki_5ht_all, drug_kaki_5ht_all);
[htpn,~,htnstats] = nansignrank(base_kaki_5ht_all, drug_kaki_5ht_all);
[~,nlp,~,nlstats] = ttest(base_kaki_nacl_all, drug_kaki_nacl_all);
[nlpn,~,nlnstats] = nansignrank(base_kaki_nacl_all, drug_kaki_nacl_all);
[~,dp,~,dstats] = ttest2(drug_kaki_5ht_all - base_kaki_5ht_all, drug_kaki_nacl_all - base_kaki_nacl_all);
[dpn,~,dnstats] = nanranksum(drug_kaki_5ht_all - base_kaki_5ht_all, drug_kaki_nacl_all - base_kaki_nacl_all);

T(:,15) = [length(find([exinfo.ismango]==0)), length(p_kaki_5ht_all), length(p_kaki_nacl_all),...
    median([base_kaki_5ht_all base_kaki_nacl_all]),...    
    median([drug_kaki_5ht_all drug_kaki_nacl_all]), htstats.tstat, htp, htnstats.signedrank, htpn, ...     
    length(find(p_kaki_5ht_all==1)),length(find(p_kaki_5ht_all==0)),...
    median(base_kaki_5ht_all), median(drug_kaki_5ht_all),...
    nlstats.tstat, nlp, nlnstats.signedrank, nlpn,  length(find(p_kaki_nacl_all==1)),length(find(p_kaki_nacl_all==0)),...
    median(base_kaki_nacl_all), median(drug_kaki_nacl_all),...
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

T(:,13) = [length(exinfo),  length(p_mango_5ht_all)+length(p_kaki_5ht_all), length(p_mango_nacl_all)+length(p_kaki_nacl_all),...
    median([base_mango_5ht_all base_mango_nacl_all base_kaki_5ht_all base_kaki_nacl_all]),...   
    median([drug_mango_5ht_all drug_mango_nacl_all drug_kaki_5ht_all drug_kaki_nacl_all]), ...
    htstats.tstat, htp, htnstats.signedrank, htpn, ...     
    length(find([p_mango_5ht_all p_kaki_5ht_all]==1)),length(find([p_mango_5ht_all p_kaki_5ht_all]==0)),...
    median([base_mango_5ht_all base_kaki_5ht_all]), median([drug_mango_5ht_all drug_kaki_5ht_all]),...
    nlstats.tstat, nlp, nlnstats.signedrank, nlpn,  ...
    length(find([p_mango_nacl_all p_kaki_nacl_all]==1)),length(find([p_mango_nacl_all p_kaki_nacl_all]==0)),...
    median([base_mango_nacl_all base_kaki_nacl_all]), median([drug_mango_nacl_all drug_kaki_nacl_all]),...
    dstats.tstat, dp, dnstats.ranksum, dpn];

if saveoption==1
    if strcmp(name, 'fixationspan')
        name = 'fixation_precision';
    end
    
    % save figure
    savefig(h, ['Z:\Corinna\SerotoninPaper\Figures\Figure08_EyeData\raw_figs\' name '.fig'])
    disp(['The figure was saved in '  ['Z:\Corinna\SerotoninPaper\Figures\Figure08_EyeData\raw_figs\' name]])
    
    % save table
    VariableNames = {'or_both','sf_both','co_both', 'sz_both',...
        'or_mango','sf_mango', 'co_mango', 'sz_mango',...
        'or_kaki','sf_kaki', 'co_kaki', 'sz_kaki','collapsed_all','mango_all','kaki_all'} ;
    RowNames = {'total', 'n_5HT', 'n_NaCl','base; median','drug; median','5HT: ttest, t = ', '5HT: ttest, p = ', '5HT: signrank, rank = ','5HT: signrank, p = ', '5ht; significant', '5ht; non-significant', '5ht_base; median','5ht; median', ...
        'nacl: ttest, t = ', 'nacl: ttest, p = ', 'nacl: signrank, rank = ','nacl: signrank, p = ', 'nacl; significant', 'nacl; non-significant',   'nacl_base; median','nacl; median',...
         'diff: ttest2, t = ', 'diff: ttest, p = ', 'diff: ranksum, rank = ','diff: ranksum, p = '};  
    T = array2table(T, 'RowNames', RowNames,'VariableNames',VariableNames);
     writetable(T, ['Z:\Corinna\SerotoninPaper\Figures\Figure08_EyeData\raw_figs\' name '.csv'],'WriteRowNames',true)
end
