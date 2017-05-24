function fitparam = fitSZ( mn, sem, sz, bootstrp )
%fitparam = fitSZ( mn, sem, sz ) 
% or
%fitparam = fitSZ( mn, sem, sz ) 
% 
% fit ratio of gaussians to size tuned spike rates
% 
% The input arguments are vectors:
% - mn:         the spike count or spike rate 
% - sem:        the standard errors of mn
% - sz:         stimulus size samples
% 
% 
% If you add a fourth argument, the tuning fit is bootstrapped and another
% field with these results is added to the output structure.
% 
% The output argument is a struct with the 
% - the fitting parameters 
% - the raw tuning curve (responses and size samples), excluding the blank:
% - well-sampled, fitted size response function 
% - the bootstrapped fitting parameters
% 
% @CL


%% parse input

noblankidx = sz<1000; % ignore blanks
sz2=sz(noblankidx); 
mn2 = mn(noblankidx); offset = min(mn2); mn2 = mn2-offset; % remove spontaneous response activity
sem2 = sem(noblankidx);


%% CL fit
[ks, kc, ws, wc, fvalcl, r2] = fitGaussRatio(sz2, mn2);


%%  assign results
fitparam.val.mn = mn(noblankidx);
fitparam.val.sem = sem2;
fitparam.val.sz = sz2;

fitparam.val.x = 0:0.03:sz2(end);
fitparam.val.y = GaussRatio(ks, kc, ws, wc, fitparam.val.x)+offset;
fitparam.wc = wc; 
fitparam.ws = ws; 
fitparam.ks = ks; 
fitparam.kc = kc; 


% the fitting quality, i.e. the explained variance
fitparam.r2 = r2;

%% compute the suppression index
% it is the response amplitude divided by the maximal response
mn4si = mn(noblankidx)-mn(~noblankidx); % remove spontaneous activity
fitparam.SI = (max(mn4si)-mn4si(end)) /  max(mn4si);

    
%% compute the preferred stimulus size
fitparam.mu = getPSZ(fitparam.val.y , fitparam.val.x, mn(noblankidx), sz2);

%% bootstrap the fitting parameters
if nargin == 4
    mn2 = mn(noblankidx);
    sem2 = sem(noblankidx);
    sz2 = sz(noblankidx);
    
    parfor i = 1:1000
        bootidx = randi(length(mn2), length(mn2), 1);
        boot(i) = fitSZ( [mn2(bootidx) mn(~noblankidx)], ...
            [sem2(bootidx) sem(~noblankidx)], ...
            [sz2(bootidx) 10001] ) ;
    end
    fitparam.boot = boot;
end

end



function prefsz = getPSZ(rfit, szfit, r, sz)
% the preferred size is the smallest radius at which
% the response exceeds 98% of the unit’s maximum response based on the
% fitted data.

[~, sortidx] = sort(sz);
r = r(sortidx);

[~, idxmax] = max(r);
if idxmax==1
    prefsz = sz(1);
else
    idx = find( rfit > 0.98*max(rfit), 1, 'first'); % get index of first response > 98%
    if isempty(idx)
        prefsz = -1;
    else
        prefsz = szfit(idx); % assign preferred size to output
    end
end

fprintf('pref sz %1.2f \n', prefsz)
end