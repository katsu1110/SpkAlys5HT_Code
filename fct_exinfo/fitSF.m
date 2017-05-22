function fitparam = fitSF( mn, sem, sf, log_flag, bootstrp)
%fitSF fits a gaussian function to the mean spike mn rates given the spatical
%frequency sf.


rng(9123234); % to always end up with the same fit in repetitive batches


% ignore blanks
i_noblank = sf<1000;
mn   = mn( i_noblank );    
sem   = sem( i_noblank );
sf  = sf( i_noblank );

if log_flag; sf = log(sf); end
    


%%% fit data to gauss
[~, sf_peak] = max(mn);
x0 = [sf_peak 4 max(mn) min(mn)]; % starting parameters (mu, sig, a, b)

if log_flag
    fo = fitoptions('Method','NonlinearLeastSquares', 'Lower', [-inf 0 0 0], ...
        'Upper', [max(sf) inf inf inf], 'StartPoint', x0); 
else
    fo = fitoptions('Method','NonlinearLeastSquares', 'Lower', [0 0 0 0], 'StartPoint', x0);
end
ft = fittype(@(mu, sig, a, b, x) gaussian(mu, sig, a, b, x), 'options', fo);
          
[f,gof2, output] = fit(sf', mn', ft);

% assign parameters
fitparam.mu = f.mu;        fitparam.sig = f.sig;
fitparam.a = f.a;          fitparam.b = f.b;
fitparam.r2 = gof2.rsquare;

% recorded data
fitparam.val.mn = mn;
fitparam.val.sem = sem;
fitparam.val.sf = sf;

% fit predicted tuning
fitparam.x = min(sf):0.01:sf(end);    
fitparam.y = gaussian(fitparam.mu, fitparam.sig, fitparam.a, fitparam.b, fitparam.x);

% % debugging
% figure; plot(fitparam.x, fitparam.y,'k'); hold on
% errorbar(sf, resp, sem, 'k');
% 

if log_flag 
    fitparam.mu = exp(fitparam.mu);
    fitparam.sig = exp(fitparam.sig);
end




if nargin == 5
    parfor i = 1:1000
        bootidx = randi(length(mn), length(mn), 1);
        boot(i) = fitSF( mn(bootidx), sem(bootidx), sf(bootidx), log_flag);
    end
    fitparam.boot = boot;
end


end


%% COST FUNCTION
function cost = cf(p, sf, resp, log_flag)
% cost function for all parameters
% note, the fit is restricted to wc<ws, i.e. the center is smaller than
% surround


mu = p(1);  sig = p(2);
a = p(3);   b = p(4);

r_pred = gaussian(mu, sig, a, b, sf);

if log_flag;     mu = exp(mu);  end


if mu<0 || sig<0 || a<0 || b<0  %restrictions
    cost = inf;
else
    if size(resp,1) ~= size(r_pred, 1)
        resp = resp';
    end
    cost = sum( (r_pred-resp).^2)  / sum( (r_pred-resp).^2); % nonlinear squared error
end
end  

