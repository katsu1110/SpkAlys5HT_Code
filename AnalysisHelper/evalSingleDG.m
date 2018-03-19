function [ex, argout] = evalSingleDG( exinfo, fname, varargin )
% batch file that calls all functions to analyse spiking activity in
% resposne to drifiting grating experiments
% 
% @CL

% accumulates the results over the course of the code
argout = {}; 

%================================================================ load data
ex = loadCluster( fname, 'ocul', exinfo.ocul, 'loadlfp', false ); % load raw data
[ex_orig, spkrate_orig] = znormex(ex, exinfo, 1);

for ps = 1:3
    ex = ex_orig;
        switch ps                
                case 1
                        suffix = '_1';
                case 2
                        suffix = '_2';
                case 3
                        suffix = '';
        end
        %============================= Analysis of stimulus selectivity using ANOVA
        if any(strcmp(varargin, 'anova')) || any(strcmp(varargin, 'all'))
            p_argout = SelectivityCheck(exinfo, ex, suffix);
            argout = [argout p_argout{:}];
        end
        %===================================================== Fano Factor analysis
        if any(strcmp(varargin, 'ff')) || any(strcmp(varargin, 'all'))
            ff_argout = FF(exinfo, ex, suffix);
            argout = [argout ff_argout{:}];
        end
        %============================================= Noise and signal correlation
        if any(strcmp(varargin, 'rsc')) || any(strcmp(varargin, 'all'))
            rsc_argout = getCoVariabilityCoeff(exinfo, fname, suffix);
            argout = [argout rsc_argout{:}];
        end
        %========================================================= Tuning curve fit
        if any(strcmp(varargin, 'tcfit')) || any(strcmp(varargin, 'all'))
            [tc_argout, tcdiff] = fitTC(exinfo, ex, fname, suffix);
            argout = [argout tc_argout{:}];
        end
        %====================================================== Tuning curve height 
        if any(strcmp(varargin, 'tcheight')) || any(strcmp(varargin, 'all'))
           argout = [argout {['tcdiff' suffix], tcdiff}];
        end
        %======================================================== Phase selectivity
        if any(strcmp(varargin, 'phasesel')) || any(strcmp(varargin, 'all'))
            switch ps
                case 1
                    ex.Trials = ex.Trials(ex.pupil_split(1).idx);
                case 2
                    ex.Trials = ex.Trials(ex.pupil_split(2).idx);
            end
            phasesel = getPhaseSelectivity(ex, 'stim', exinfo.param1);
            argout = [argout {['phasesel' suffix], phasesel}];
        end
        % ==================================== fixation precision & microsaccade
        if any(strcmp(varargin, 'phasesel')) || any(strcmp(varargin, 'all'))
                fixspan = fixationSpan(ex, 0);
                amp = [];
                peakv = [];
                duration = [];
                angle = [];
                for i = 1:length(fixspan.trial)
                    amp = [amp fixspan.trial(i).microsaccade.amp];
                    peakv = [peakv fixspan.trial(i).microsaccade.peakv];
                    duration = [duration fixspan.trial(i).microsaccade.duration];
                    angle = [angle fixspan.trial(i).microsaccade.angle];
                end
                argout = [argout {['fixationspan' suffix],  [fixspan.accuracy, fixspan.theta, fixspan.SI]', ...
                        ['variance_eye' suffix], [fixspan.variance]',...
                        ['recentering_win' suffix],  [fixspan.recentering.window_offsetX; fixspan.recentering.window_offsetY]', ...
                        ['recentering_eye' suffix], [fixspan.recentering.eye_offsetX; fixspan.recentering.eye_offsetY]', ...
                        ['microsac_amplitude' suffix], amp, ...
                        ['microsac_peakv' suffix], peakv, ...
                        ['microsac_duration' suffix], duration, ...
                        ['microsac_angle' suffix], angle, ...
                        ['microsac_counts' suffix], arrayfun(@(x) x.microsaccade.counts, fixspan.trial)}];
        end
        %===================================== organize results as output arguments
        ex = ex_orig;
        if ismember(ps, [1 2])     
             switch ps
                    case 1
                        ex.Trials = ex.Trials(ex.pupil_split(1).idx);
                    case 2
                        ex.Trials = ex.Trials(ex.pupil_split(2).idx);
             end
            [ex, spkrate] = znormex(ex, exinfo, 1);
        else
            spkrate = spkrate_orig;
        end
        argout =  [argout {['rateMN' suffix], [spkrate.mn]', ['rateVARS' suffix], [spkrate.var]', ...
            ['ratePAR' suffix], [spkrate.(exinfo.param1)]', ['rateSME' suffix], [spkrate.sem]', ...
            ['rateSD' suffix], [spkrate.sd]', ['rawspkrates' suffix], [spkrate.raw]', ...
            ['rate_resmpl' suffix], {spkrate.resamples}, ['nrep' suffix], [spkrate.nrep]',...
             }];
end

end
    
    
%% subfunctions to compute the analysis

%% co-variability

function res_argout = getCoVariabilityCoeff(exinfo, fname, suffix)
% compute the noise and signal correlation between single and multiunit
% activity based on the full and second half of the experiment.
% Reducing the analysis window to the second half reduces the risk of
% artificially increased co-variability by effects of adaptation or the
% gradual increase of the 5HT effect in the beginning.
% 


% load data
ex_su = loadCluster(fname, 'ocul', exinfo.ocul, 'loadlfp', false ); % cluster 1 or 2 <=> single unit activity
[ex_su] = znormex(ex_su, exinfo, 1);

fname_mu = strrep(fname, exinfo.cluster, 'c0'); % cluster 0 <=> multiunit activity
ex_mu = loadCluster( fname_mu, 'ocul', exinfo.ocul, 'loadlfp', false );

switch suffix
    case '_1'
        ex_su.Trials = ex_su.Trials(ex_su.pupil_split(1).idx);
        ex_mu.Trials = ex_mu.Trials(ex_su.pupil_split(1).idx);
    case '_2'
        ex_su.Trials = ex_su.Trials(ex_su.pupil_split(2).idx);
        ex_mu.Trials = ex_mu.Trials(ex_su.pupil_split(2).idx);
end

%%% all trials
[rsc, prsc, rsig, prsig, c0geomn, trials_su, trials_mu] = ...
    getCoVariabilityCoeff_helper(ex_su, ex_mu, exinfo);
res_fullexp = {['rsc' suffix], rsc, ['prsc' suffix], prsc, ['rsig' suffix], rsig,['prsig' suffix], prsig, ...
    ['c0geomn' suffix], c0geomn, ['trials_c0' suffix], trials_mu, ['trials_c1' suffix], trials_su};


if isempty(suffix)
    %%% 2nd half
    ex_su.Trials = getPartialTrials(ex_su.Trials, 2); 
    ex_mu.Trials = getPartialTrials(ex_mu.Trials, 2);

    [rsc_2nd, prsc_2nd, rsig_2nd, prsig_2nd, c0geomn_2nd] = ...
        getCoVariabilityCoeff_helper(ex_su, ex_mu, exinfo);

    res_2ndhalf = {['rsc_2nd' suffix], rsc_2nd, ['prsc_2nd' suffix], prsc_2nd, ['rsig_2nd' suffix], rsig_2nd, ...
        ['prsig_2nd' suffix], prsig_2nd, ['c0geomn_2nd' suffix], c0geomn_2nd};
else
    res_2ndhalf = [];
end

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

[ex_su, spkrate_su] = znormex(ex_su, exinfo, 0);
[ex_mu, spkrate_mu] = znormex(ex_mu, exinfo, 0); % z norm

% Noise Correlation
[rsc, prsc] = corr([ex_su.Trials.zspkcount]', [ex_mu.Trials.zspkcount]', ...
    'rows', 'complete');

% Signale Correlation
[rsig, prsig] = corr([spkrate_su.cmn]', [spkrate_mu.cmn]', ...
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
function argout = SelectivityCheck(exinfo, ex, suffix)
% test whether the unit shows significant selectivity to a stimulus by
% performing an ANOVA on the averaged spike rates of the different stimulus
% conditions
%
% ignore trials with blanks for this

switch suffix
    case '_1'
        ex.Trials = ex.Trials(ex.pupil_split(1).idx);
    case '_2'
        ex.Trials = ex.Trials(ex.pupil_split(2).idx);
end

idx_nonblank = getNoBlankIdx(exinfo, ex);

p_anova = anova1([ex.Trials(idx_nonblank).spkRate], ...
    [ex.Trials(idx_nonblank).(exinfo.param1)],'off');

argout  = {['p_anova' suffix], p_anova};

end


%%
function argout = FF(exinfo, ex, suffix)
% computation of spiking variability 
% this function calls FanoFactors.mat, a function that computes different
% metrices of variability (classic FF as variance/mean, binned FF, etc.)
% 
% as for the co-variability the analysis is once based on the entire
% experiment and once based on the second half only.
% 
% 
%

switch suffix
    case '_1'
        ex.Trials = ex.Trials(ex.pupil_split(1).idx);
    case '_2'
        ex.Trials = ex.Trials(ex.pupil_split(2).idx);
end

% results for the full experiment
[ex, spkrate] = znormex(ex, exinfo, 0);

[ ff.classic, ff.fit, ff.mitchel, ff.church ] = ...
    FanoFactors( ex, [spkrate.cmn], [spkrate.cvar], [spkrate.nrep], exinfo.param1);

if isempty(suffix)
    % results for the second experiment
    ex.Trials = getPartialTrials(ex.Trials, 2);
    [ex, spkrate] = znormex(ex, exinfo, 0);

    [ ff.classic_2ndhalf, ff.fit_2ndhalf, ff.mitchel_2ndhalf, ff.church_2ndhalf ] = ...
        FanoFactors( ex, [spkrate.cmn], [spkrate.cvar], [spkrate.nrep], exinfo.param1);
end

argout = {['ff' suffix], ff};
end


%%
function [argout, tcdiff] = fitTC(exinfo, ex, fname, suffix)
% fitting descriptive tuning functions to the raw tuning curves
%
%
% each feature (orientation, contrast, ...) has it's descriptive function
% (e.g. gaussian, sigmoid, ...). We therefore have to distinguish between
% features and call different fitting functions.
% 

switch suffix
    case '_1'
        ex.Trials = ex.Trials(ex.pupil_split(1).idx);
    case '_2'
        ex.Trials = ex.Trials(ex.pupil_split(2).idx);
end

[ex, spkrate] = znormex(ex, exinfo, 0);

% tuning curve height
minspk = min([spkrate.mn]);
maxspk = max([spkrate.mn]);
tcdiff = (maxspk - minspk) / mean([maxspk, minspk]);

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

argout = {['fitparam' suffix], fitparam};

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