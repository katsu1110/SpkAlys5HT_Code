function [ex, argout] = evalSingleDG( exinfo, fname, varargin )
% batch file that calls all functions to analyse spiking activity in
% resposne to drifiting grating experiments
% 
% @CL

% accumulates the results over the course of the code
argout = {}; 

%================================================================ load data
ex = loadCluster( fname, 'ocul', exinfo.ocul, 'loadlfp', false ); % load raw data
[ex, spkrate] = znormex(ex, exinfo);

%==================================================== Direction Selectivity
if any(strcmp(varargin, 'dirsel')) || any(strcmp(varargin, 'all'))
    
    if strcmp(exinfo.param1, 'or')
        dirsel = getDIRsel(fname);
    else
        dirsel = -1;
    end
    
    argout = [argout {'dirsel', dirsel}];
end

%======================================================== Circualr Variance
if any(strcmp(varargin, 'circvar')) || any(strcmp(varargin, 'all'))
    
    if strcmp(exinfo.param1, 'or')
        i_theta = [spkrate.or]<360;
        circvar = compCircularVariance([spkrate(i_theta).mn], [spkrate(i_theta).or] );
    else
        circvar = -1;
    end
    
    argout = [argout {'circvar', circvar}];
end


%============================= Analysis of stimulus selectivity using ANOVA
if any(strcmp(varargin, 'anova')) || any(strcmp(varargin, 'all'))
    p_argout = SelectivityCheck(exinfo, ex);
    argout = [argout p_argout{:}];
end
%===================================================== Fano Factor analysis
if any(strcmp(varargin, 'ff')) || any(strcmp(varargin, 'all'))
    ff_argout = FF(exinfo, ex);
    argout = [argout ff_argout{:}];
end
%============================================= Noise and signal correlation
if any(strcmp(varargin, 'rsc')) || any(strcmp(varargin, 'all'))
    rsc_argout = getCoVariabilityCoeff(exinfo, fname);
    argout = [argout rsc_argout{:}];
end
%========================================================= Tuning curve fit
if any(strcmp(varargin, 'tcfit')) %|| any(strcmp(varargin, 'all'))
    tc_argout = fitTC(exinfo, spkrate, ex, fname);
    argout = [argout tc_argout{:}];
end
%====================================================== Tuning curve height 
if any(strcmp(varargin, 'tcheight')) || any(strcmp(varargin, 'all'))
   minspk = min([spkrate.mn]);
   maxspk = max([spkrate.mn]);
   tcdiff = (maxspk - minspk) / mean([maxspk, minspk]);
   argout = [argout {'tcdiff', tcdiff}];
end
%======================================================== Phase selectivity
if any(strcmp(varargin, 'phasesel')) || any(strcmp(varargin, 'all'))
    phasesel = getPhaseSelectivity(ex, 'stim', exinfo.param1);
    argout = [argout {'phasesel', phasesel}];
end
%===================================== organize results as output arguments
argout =  [argout {'rateMN', [spkrate.mn]', 'rateVARS', [spkrate.var]', ...
    'ratePAR', [spkrate.(exinfo.param1)]', 'rateSME', [spkrate.sem]', ...
    'rateSD', [spkrate.sd]', 'rawspkrates', [spkrate.sd]', ...
    'rate_resmpl', {spkrate.resamples}, ...
    'nrep', [spkrate.nrep]'}];

end
    
    
%% subfunctions to compute the analysis

%% co-variability

function res_argout = getCoVariabilityCoeff(exinfo, fname)
% compute the noise and signal correlation between single and multiunit
% activity based on the full and second half of the experiment.
% Reducing the analysis window to the second half reduces the risk of
% artificially increased co-variability by effects of adaptation or the
% gradual increase of the 5HT effect in the beginning.
% 

% load data
ex_su = loadCluster(fname, 'ocul', exinfo.ocul, 'loadlfp', false ); % cluster 1 or 2 <=> single unit activity

fname_mu = strrep(fname, exinfo.cluster, 'c0'); % cluster 0 <=> multiunit activity
ex_mu = loadCluster( fname_mu, 'ocul', exinfo.ocul, 'loadlfp', false );


%%% all trials
[rsc, prsc, rsig, prsig, c0geomn, trials_su, trials_mu] = ...
    getCoVariabilityCoeff_helper(ex_su, ex_mu, exinfo);
res_fullexp = {'rsc', rsc, 'prsc', prsc, 'rsig', rsig,'prsig', prsig, ...
    'c0geomn', c0geomn, 'trials_c0', trials_mu, 'trials_c1', trials_su};


%%% 2nd half
ex_su.Trials = getPartialTrials(ex_su.Trials, 2); 
ex_mu.Trials = getPartialTrials(ex_mu.Trials, 2);

[rsc_2nd, prsc_2nd, rsig_2nd, prsig_2nd, c0geomn_2nd] = ...
    getCoVariabilityCoeff_helper(ex_su, ex_mu, exinfo);

res_2ndhalf = {'rsc_2nd', rsc_2nd, 'prsc_2nd', prsc_2nd, 'rsig_2nd', rsig_2nd, ...
    'prsig_2nd', prsig_2nd, 'c0geomn_2nd', c0geomn_2nd};


% concatenate results
res_argout = [res_2ndhalf, res_fullexp];
    
end


function [rsc, prsc, rsig, prsig, c0geomn, trials_su, trials_mu] = ...
    getCoVariabilityCoeff_helper(ex_su, ex_mu, exinfo)
% compute the noise correlation and signal correlation
%
% noise correlations are computed as the pearson correlation between
% stimulus z-normed spike counts.
%
% signal correlations are computed as the pearson correlation between raw
% tuning curves.

[ex_su, spkrate_su ,spkcount_su] = znormex(ex_su, exinfo);
[ex_mu, spkrate_mu, spkcount_mu] = znormex(ex_mu, exinfo); % z norm

% Noise Correlation
[rsc, prsc] = corr([ex_su.Trials.zspkcount]', [ex_mu.Trials.zspkcount]', ...
    'rows', 'complete');

% Signale Correlation
[rsig, prsig] = corr([spkcount_su.mn]', [spkcount_mu.mn]', ...
    'rows', 'complete');


trials_su = struct('spkRate', [ex_su.Trials.spkRate], ...
    'zspkCount', [ex_su.Trials.zspkcount], ...
    'param',[ex_su.Trials.(exinfo.param1)]);
trials_mu= struct('spkRate', [ex_mu.Trials.spkRate], ...
    'zspkCount', [ex_mu.Trials.zspkcount], ...
    'param',[ex_mu.Trials.(exinfo.param1)]);

% the geometric mean of spike rates
c0geomn = geomean([ mean([spkrate_mu.mn]), mean([spkrate_su.mn]) ]);

end




%%
function argout = SelectivityCheck(exinfo, ex)
% test whether the unit shows significant selectivity to a stimulus by
% performing an ANOVA on the averaged spike rates of the different stimulus
% conditions
%
% ignore trials with blanks for this

idx_nonblank = getNoBlankIdx(exinfo, ex);

p_anova = anova1([ex.Trials(idx_nonblank).spkRate], ...
    [ex.Trials(idx_nonblank).(exinfo.param1)],'off');

argout  = {'p_anova', p_anova};

end


%%
function argout = FF(exinfo, ex)
% computation of spiking variability 
% this function calls FanoFactors.mat, a function that computes different
% metrices of variability (classic FF as variance/mean, binned FF, etc.)
% 
% as for the co-variability the analysis is once based on the entire
% experiment and once based on the second half only.
% 
% 
%


% results for the full experiment
[ex, spkrate, spkcount] = znormex(ex, exinfo);

[ ff.classic, ff.fit, ~, ~ ] = ...
    FanoFactors( ex, [spkcount.(exinfo.param1)], ...
    [spkcount.mn], [spkcount.var], [spkrate.nrep], exinfo.param1);


% results for the experiment neglecting the first 20 trials
ex_20plus = ex;
ex_20plus.Trials = ex_20plus.Trials(21:end);

[ex_20plus, spkrate_20plus, spkcount_20plus] = znormex(ex_20plus, exinfo);

[ ff.classic_20plus, ff.fit_20plus, ~, ~ ] = ...
    FanoFactors( ex_20plus,  [spkcount_20plus.(exinfo.param1)], ...
    [spkcount_20plus.mn], [spkcount_20plus.var], [spkrate_20plus.nrep], exinfo.param1);


% results for the 2nd half of the experiment
ex_2ndhalf = ex;
ex_2ndhalf.Trials = getPartialTrials(ex_2ndhalf.Trials, 2);
[ex_2ndhalf, spkrate_2ndhalf, spkcount_2ndhalf] = znormex(ex_2ndhalf, exinfo);

[ ff.classic_2ndhalf, ff.fit_2ndhalf, ~, ~ ] = ...
    FanoFactors( ex_2ndhalf,  [spkcount_2ndhalf.(exinfo.param1)], ...
    [spkcount_2ndhalf.mn], [spkcount_2ndhalf.var], [spkrate_2ndhalf.nrep], exinfo.param1);

argout = {'ff', ff};
end

%%
function argout = fitTC(exinfo, spkrate, ex, fname)
% fitting descriptive tuning functions to the raw tuning curves
%
%
% each feature (orientation, contrast, ...) has it's descriptive function
% (e.g. gaussian, sigmoid, ...). We therefore have to distinguish between
% features and call different fitting functions.
% 

if strcmp(exinfo.param1, 'or')
    % fitting the orientation tuning curve
    
    fitparam = fitOR( [spkrate.mn], [spkrate.sem], [spkrate.or] );
    
elseif strcmp(exinfo.param1, 'co')
    % fitting the contrast tuning curve
    
    fitparam = fitCO([spkrate.mn], [spkrate.(exinfo.param1)]);
    
    % check for underspampled contrast variations
    % in many experiments the area around the putative c50 is sparsely
    % sampled. in order to make a useful statement about c50 or saturation,
    % we need to check, whether we also have selectivity before and after
    % c50. I do so by computing the ANOVA on the response to contrasts
    % smaller or greater than the putative c50. 
    idx_nonblank = getNoBlankIdx(exinfo, ex); % ignore the blank trials

    i_cosmc50 = [ex.Trials.(exinfo.param1)]<=fitparam.c50 & idx_nonblank; % smaller than c50
    fitparam.undersmpl(1) = anova1([ex.Trials(i_cosmc50).spkRate], [ex.Trials(i_cosmc50).(exinfo.param1)],'off');
    
    i_cohic50 = [ex.Trials.(exinfo.param1)]>fitparam.c50 & idx_nonblank;% higher than c50
    fitparam.undersmpl(2) = anova1([ex.Trials(i_cohic50).spkRate], [ex.Trials(i_cohic50).(exinfo.param1)],'off');
    
    
    % in case of the drug file, test the contrast gain and activity gain model 
    if strfind( fname , exinfo.drugname )
        fitparam.sub = fitCO_drug([spkrate.mn], [spkrate.(exinfo.param1)], exinfo.fitparam);
    end
    
elseif strcmp(exinfo.param1, 'sz')
    % fitting the size tuning curve
    fitparam = fitSZ([spkrate.mn], [spkrate.sem], [spkrate.(exinfo.param1)]);
    
elseif strcmp(exinfo.param1, 'sf')
    % fitting the spatial frequency tuning curve
    fitparam_lin = fitSF([spkrate.mn], [spkrate.sem], [spkrate.(exinfo.param1)], false );
    fitparam_log = fitSF([spkrate.mn], [spkrate.sem], [spkrate.(exinfo.param1)], true );
    
    % first entry linear fit, second entry log scaled fit
    fitparam.others = {fitparam_lin; fitparam_log};
else
    fitparam = [];
end

argout = {'fitparam', fitparam};

end

%%
function idx_nonblank = getNoBlankIdx(exinfo, ex)
% trials with blanks

if strcmp(exinfo.param1, 'co') %sanity check, in case the blank was coded as 0% contrast
    idx_nonblank = [ex.Trials.(exinfo.param1)] ~= ex.exp.e1.blank & ...
        [ex.Trials.(exinfo.param1)] ~= 0;
else
    idx_nonblank = [ex.Trials.(exinfo.param1)] < ex.exp.e1.blank;
end
end