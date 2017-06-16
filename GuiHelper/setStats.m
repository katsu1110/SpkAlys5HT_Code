function stats = setStats( h, inclusioncrit)
% adds a string to the figure handle h's UserData with all statistical
% tests and their results for the attached data set. The string can be
% called via h.UserData.Stats
% 
% It contains the two-sample results (5HT vs. NaCl) and the paired test
% results for each condition (Baselin vs 5HT, Baseline vs NaCl) for both
% and each individual monkey.
%

dat = h.UserData;
idx = [dat.expInfo.is5HT];

x = dat.x; y = dat.y;

%%% all
s = getStatsHelper(x, y, idx);

%%% mango
idx2 =[dat.expInfo.ismango];
s_ma = getStatsHelper(x(idx2), y(idx2), idx(idx2));

%%% kaki
idx2 =~[dat.expInfo.ismango];
s_ka = getStatsHelper(x(idx2), y(idx2), idx(idx2));



% concatenate and save in UserData
dat.Stats = sprintf(['X = ' dat.xlab ' \nY = ' dat.ylab '\n\n' inclusioncrit '\n\nALL' s...
    '\n----------\n\n' 'MANGO' s_ma ...
    '\n----------\n\n' 'KAKI' s_ka]);
set(h, 'UserData', dat);
stats = dat.Stats;

end


function s = getStatsHelper(X, Y, idx)
% calls the statistics for a particular data set
% idx is the index indicating 5HT (==1) and NaCl (==0) experiments

s_5HT = getStatsPairedSample(X(idx), Y(idx));
s_5HT_corr = getCorr(X(idx)', Y(idx)');
s_NaCl = getStatsPairedSample(X(~idx), Y(~idx));
s_NaCl_corr = getCorr(X(~idx)', Y(~idx)');

s_x = getStatsTwoSample(X(idx), X(~idx));
s_y = getStatsTwoSample(Y(idx), Y(~idx));
s_diff = getStatsTwoSample(X(idx)-Y(idx), X(~idx)-Y(~idx));


s = ['\n----5HT\n' s_5HT s_5HT_corr...
    '\n----NaCl\n' s_NaCl s_NaCl_corr...
    '\n----5HT vs. NaCl\n x: ' s_x  ' y: ' s_y  ' x-y: ' s_diff];
end


function s = getCorr(X, Y)
% correlation between x and y, both parametric (Pearson) and non-parametric
% (Spearman)
if isempty(X)
    rho_p=0; p_p=1;
    rho_sp=0; p_sp=1;
else
    [rho_p, p_p] = corr(X, Y, 'type', 'Pearson');
    [rho_sp, p_sp] = corr(X, Y, 'type', 'Spearman');
end
s = sprintf('correlation Pearson: rho=%1.2f  p=%1.2e, \t Spearman: rho=%1.2f  p=%1.2e \n', ...
        rho_p, p_p, rho_sp, p_sp);
    
end

function s = getStatsTwoSample(A, B)
% two-sample tests parametric and non-parametric
if isempty(A) || isempty(B)
    ptt = 1;
    pwil = 1;
else
    % two sample comparison
    [~, ptt] = ttest2(A, B);
    pwil = ranksum(A, B);
end

s = sprintf('2-sample ttest p=%1.2e, \t wilcoxon p=%1.2e \n', ptt, pwil);

end

function s = getStatsPairedSample(X, Y)
% paired test parametric and non-parametric

s_base  = getDistParam(X);
s_drug = getDistParam(Y);
s_diff = getDistParam(X-Y);

% check for normal distribution
if isempty(X)
    h = 0;
    ptt = 1;
    psr = 1;
else
    h = kstest(X-Y);
    [~, ptt] = ttest(X, Y);
    psr = signrank(X, Y);
end

s = sprintf(['X' s_base 'Y' s_drug 'Diff X-Y' s_diff ...
    'normal dist X-Y h=%1.0f, paired t-test p=%1.2e, \t paired signrank p=%1.2e \n\n' ], ...
    h, ptt, psr);

end



function s = getDistParam(A)
% descriptive and normative statistics of a data set A

n = length(A);

% distribution values
if length(A)>2
    b = bootstrp(1000,@median,A);
    prct =  prctile(b, [5 50 95]);
else
    prct = nan(3,1);
end
prctl(2) = median(A);

mn_ = nanmean(A);
std_ = nanstd(A);

if any(A<0)
    geomn_ = -inf;
else
    geomn_ = geomean(A);
end

% check for normal distribution
if isempty(A)
    psphericity = nan;
    pttest = nan;
    psignr = nan;
else 
    [~,psphericity] = kstest(A);
    
    % test for significant difference to normal distribution 
    [~, pttest] = ttest(A, 0);   
    psignr = signrank(A);
end




s = sprintf(['(N = %1.0f) \n' ...
    'percentile (5, 50, 95): %1.2e/%1.2e/%1.2e, \t mean: %1.2e +- %1.2f SD, \t GM: %1.2e \n' ...
    'test for normal dist p=%1.2e, \t t-test vs. 0 p = %1.2e, \t signrank test vs 0 p: %1.2e \n\n'], ...
    n, prct, mn_, std_, geomn_, psphericity, pttest, psignr);

end
