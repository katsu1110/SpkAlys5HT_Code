function fitparam = fitSZ( mn, sem, sz, bootstrp )
% assigns the critical parameters and estimates to the fitparam structure
%
%


idx = sz<1000;
sz2=sz(idx); 
mn2 = mn(idx); offset = min(mn2); mn2 = mn2-offset;
sem2 = sem(idx);

sz_4fit = sz2;


% CL fit
[ks, kc, ws, wc, fvalcl, r2] = fitGaussRatio(sz_4fit, mn2);


%  assign results
fitparam.val.mn = mn(idx);
fitparam.val.sem = sem2;
fitparam.val.sz = sz2;

fitparam.val.x = 0:0.03:sz2(end);
fitparam.val.y = GaussRatio(ks, kc, ws, wc, fitparam.val.x)+offset;
fitparam.wc = wc; 
fitparam.ws = ws; 
fitparam.ks = ks; 
fitparam.kc = kc; 

mn4si = mn(idx)-mn(~idx);
fitparam.SI = (max(mn4si)-mn4si(end)) /  max(mn4si);

    

fitparam.mu = getPSZ(fitparam.val.y , fitparam.val.x, mn(idx), sz2);
fitparam.r2 = r2;

if nargin == 4
    mn2 = mn(idx);
    sem2 = sem(idx);
    sz2 = sz(idx);
    
    parfor i = 1:1000
        bootidx = randi(length(mn2), length(mn2), 1);
        boot(i) = fitSZ( [mn2(bootidx) mn(~idx)], ...
            [sem2(bootidx) sem(~idx)], ...
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


% 
%  
% function Q = MyRatioOfGaussians(X0,x)
% 
% f1 = @(t1) exp(-(t1/X0(1)).^2);
% f2 = @(t2) exp(-(t2/X0(2)).^2);
% for n=1:length(x)
% Q(n) = X0(3)*(2/sqrt(pi)*integral(f1,0,x(n)))^2/...
%     (1+(2/sqrt(pi)*X0(4)*integral(f2,0,x(n)))^2);
% end
% 
% end
%  
 