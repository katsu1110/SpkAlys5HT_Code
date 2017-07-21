function ds = getDIRsel( exp_fname, varargin)
% ds = getDIRsel( exp_fname, varargin )
%
% computes the direction selectivity from the experiment.
% exp_fname is the filename of the corresponding ex file.
% 
% ds is the direction selectivity computed as (Rpref-Rnull)/(Rpref+Rnull)
% with Rpref being the maximum response (at orientation theta_pref) and
% Rnull the response to theta_pref+180.
%
%
% see Prince et al. (2000) or Ruff and Cohen (2016)
% http://jn.physiology.org/content/jn/87/1/191.full.pdf
% http://www.jneurosci.org/content/jneuro/36/28/7523.full.pdf
% 
%@CL, 29.06.2017

if any(strcmp(varargin, 'plot'))
    p_flag = 1;
else 
    p_flag = 0;
end


% load the experiment's recording data
load(exp_fname);
addspkRate;

% use only rewarded trials
Trials = ex.Trials([ex.Trials.Reward]==1); 

% compute the spontaneouse activity
Rblank = mean ([Trials([Trials.st]==1).spkRate]);

% stimulus samples (i.e. the stimulus orientations)
theta = unique([Trials([Trials.st]==1).or]);


% average across stimulus conditions
for n = 1:length(theta)
    R_raw(n) = mean( [Trials([Trials.or]==theta(n)).spkRate] );
end

R = R_raw;

%R = R_raw-Rblank;
% R(R<0) = 0; % in case the spontaneous activity is higher than driven responses

% response to the preferred direction
[~, maxi] = max(R);
Rpref = R(maxi);


% response to the null direction (that is the reverse drifting direction)
theta_null = mod(theta(maxi)+180, 360);

Rnull = R( theta == theta_null );

% finally, the direction selectivity is 1 minus the ratio between the
% response difference and the preferred direction response
ds = (Rpref-Rnull)/Rpref;
 

if isempty(Rnull)
    ds = nan;
    warning('there is no null direction. direction selecitivy was set to nan');
end



if p_flag
    
   plot(theta, R_raw, '--'); hold on;
   plot([theta(1), theta(end)], ones(1,2)*Rblank, '--');
    
   scatter(theta(maxi), max(R_raw), 'k', 'filled', 'Marker', 'v');
   scatter(theta_null, max(R_raw),  'k', 'Marker', 'v');
end
    
end

