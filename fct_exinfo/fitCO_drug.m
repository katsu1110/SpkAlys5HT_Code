function fitparam =  fitCO_drug(resp, co, basefit)
% finds the best fit for the 

co(co>1000) =0;
ss_tot = sum( (resp - mean(resp)).^2 );
fitparam = basefit;


a0 = 1;
fitparam.a_cg       = fminsearch(@(a) contrastgain(co, resp, basefit, a), a0);
ss_res_cg           = contrastgain(co, resp, basefit, fitparam.a_cg);
fitparam.r2_cg      = 1 - (ss_res_cg/ss_tot); 


fitparam.a_ag       = fminsearch(@(a) activitygain(co, resp, basefit, a), a0);
ss_res_ag           = activitygain(co, resp, basefit, fitparam.a_ag);
fitparam.r2_ag      = 1 - (ss_res_ag/ss_tot); 

fitparam.a_rg       = fminsearch(@(a) responsegain(co, resp, basefit, a), a0);
ss_res_rg           = responsegain(co, resp, basefit, fitparam.a_ag);
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





function ss = contrastgain(co, r, basefit, a)
% activity gain function, returns sum of squares of residuals
% a is the modulation factor modulating c50
    
n = basefit.n;
rpred = basefit.rmax .* ( co.^n ./ (co.^n + (a* (basefit.c50^n)) ) ) + basefit.m;
    
ss = sum((rpred-r).^2);

if a<0;  ss = 10^6;  end

end



function ss = activitygain(co, r, basefit, a)
% activity gain function, returns sum of squares of residuals
% a is the modulation factor of the entire function
    
n = basefit.n;
rpred = a .* ( basefit.rmax .* ( co.^n ./ (co.^n + (basefit.c50^n) ) ) + basefit.m);
    
ss = sum((rpred-r).^2);
if a<0;  ss = 10^6;  end
end



function ss = responsegain(co, r, basefit, a)
% activity gain function, returns sum of squares of residuals
% a is the modulation factor of the entire function except the spontanouse
% firing

n = basefit.n;
rpred = a .* ( basefit.rmax .* ( co.^n ./ (co.^n + (basefit.c50^n) ) ) )+ basefit.m;
    
ss = sum((rpred-r).^2);
if a<0;  ss = 10^6;  end
end