function fitparam =  fitCO_drug(mn, co, fitparam)
% fitparam =  fitCO_drug(resp, co, fitparam)
% 
% test different modulations of the fitted baseline tuning curve (contrast
% gain model, activity gain model, response gain model) and extent the
% fitparam structure witht the results.
% 
% input:
% mn        - mean stimulus response in the drug conditio
% co        - contrast samples corresponding to mn
% fitparam  - results of the tuning fit fo the baseline experiment ( see
%               fitCO for more information)
% 
% 
% outpu: 
% fitparam with the additiona fields:
% .a_cg, .r2_ag     - mbest fitting modulation factor in the contrast gain
%       model and the correspongin fitting quality (i.e. explained variance)
% 
% .a_ag, .r2_ag     - same for the acticity gain model
% 
% .a_cg, .r2_ag     - same for the response gain model
% 
% see Williford and Maunsell, 2006
% 
% 
% @CL


idxblank = co>1 | co==0;
co = co(~idxblank); % safety check, ignore blanks
mn = mn(~idxblank);

 
% in all cases, we start with a model that has no change
a0 = 1; %starting value. 
ss_tot = sum( (mn - mean(mn)).^2 ); % sum of squared errors 


% contrast gain model, scaling c50
fitparam.a_cg       = fminsearch(@(a) contrastgain(co, mn, fitparam, a), a0);
ss_res_cg           = contrastgain(co, mn, fitparam, fitparam.a_cg);
fitparam.r2_cg      = 1 - (ss_res_cg/ss_tot); 

% activity gain model, scaling the enitre function
fitparam.a_ag       = fminsearch(@(a) activitygain(co, mn, fitparam, a), a0);
ss_res_ag           = activitygain(co, mn, fitparam, fitparam.a_ag);
fitparam.r2_ag      = 1 - (ss_res_ag/ss_tot); 

% activity gain model, scaling the response gain
fitparam.a_rg       = fminsearch(@(a) responsegain(co, mn, fitparam, a), a0);
ss_res_rg           = responsegain(co, mn, fitparam, fitparam.a_ag);
fitparam.r2_rg      = 1 - (ss_res_rg/ss_tot); 



% find the best fit and reset the fitting parameters accordingly
r2 = [fitparam.r2_cg, fitparam.r2_ag, fitparam.r2_rg];
[~, maxi] = max(r2);

if maxi ==1 
    fitparam.c50 = fitparam.c50 * fitparam.a_cg;    
elseif maxi == 2
    fitparam.rmax = fitparam.rmax * fitparam.a_ag;
    fitparam.m = fitparam.m * fitparam.a_ag;
elseif maxi == 3
    fitparam.rmax = fitparam.rmax * fitparam.a_rg;
end

fitparam.r2 = r2(maxi);

fitparam.x = eps:0.001:1;
fitparam.y = hyperratiofct(fitparam.x, fitparam.rmax, fitparam.c50, fitparam.n, fitparam.m);
fitparam.auc = sum(hyperratiofct( 0.1:0.1:1, fitparam.rmax, fitparam.c50, fitparam.n, fitparam.m));
end





function ss = contrastgain(co, r, fitparam, a)
% contrast gain function, returns sum of squared residuals
% a is the modulation factor modulating c50
    
n = fitparam.n;
rpred = fitparam.rmax .* ( co.^n ./ (co.^n + (a* (fitparam.c50^n)) ) ) + fitparam.m;
    
ss = sum((rpred-r).^2);

if a<0;  ss = 10^6;  end

end



function ss = activitygain(co, r, fitparam, a)
% activity gain function, returns sum of squared residuals
% a is the modulation factor of the entire function
    
n = fitparam.n;
rpred = a .* ( fitparam.rmax .* ( co.^n ./ (co.^n + (fitparam.c50^n) ) ) + fitparam.m);
    
ss = sum((rpred-r).^2);
if a<0;  ss = 10^6;  end
end



function ss = responsegain(co, r, fitparam, a)
% response gain function, returns sum of squared residuals
% a is the modulation factor of the entire function except the spontanouse
% firing

n = fitparam.n;
rpred = a .* ( fitparam.rmax .* ( co.^n ./ (co.^n + (fitparam.c50^n) ) ) )+ fitparam.m;
    
ss = sum((rpred-r).^2);
if a<0;  ss = 10^6;  end
end