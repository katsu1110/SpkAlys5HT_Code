function tcfit = fitCO(mn, co_all, bootstrp)
% fitparam = fitCO(mn, co_in)
% or
% fitparam = fitCO(mn, co_in, bootstrp)
% 
% fit hyperbolic function to contrast spike rates
% 
% The input arguments are vectors:
% - mn:         the spike count or spike rate 
% - co_all:     stimulus contrast samples
% 
% 
% The fit is restricted in the subfunciton cf().
% If you add a fourth argument, the tuning fit is bootstrapped and another
% field with these results is added to the output structure.
%
% 
% The output argument is a struct with the 
%  - fitting results: 
%       rmax (max. response), m (spontaneous resp), n (steepness), c50
%       (contrast at inflection point), the of explained variance of the
%       fit r2
% - raw tuning curve, excluding the blank:
%       val (mean firing rate), co (stimulus samples)
% - well-sampled, fitted contrast response function 
%       x (in units of mn) and y (in degree)
%
%
% 
% @CL 21.06.2016


%%% the random number generator seed
%%% do not change 
rng(9123234); 


%%% neglect spontaneous activity in the input arguments
co = co_all(co_all<=1);
mn = mn(co_all<=1);


%%% fitting settings
opt = optimset('MaxIter', 10^4, 'MaxFunEvals', 10^4, 'TolFun', 0.001,'Display','off');


%%% in order to improve the fitting algorithm, I repeat the procedure with
%%% varying starting points
n_start = 20; 
p0 = [randn(n_start, 1).* repmat(max(mn), n_start, 1), ...  % rmax 
    abs(randn(n_start,1))+0.1, ...                    % c50
    abs(randn(n_start,1)), ...                     % n
    abs(randn(n_start,1))*5+min(mn)];         % m
p0(p0(:,2)<0, 2) = 0; % don't start with c50 smaller than 0
p0(p0(:,2)>1, 2) = 1; % don't start with c50>1

% fitting
parfor i =1:n_start
    param(:, i)= fminsearch(@(p) cf(co, mn, p), p0(i,:), opt);
    ss(i) = cf(co, mn, param(:, i));
end

%%% choose the best fit according to the squared error
[~, mini] = min(ss);
ss = ss(mini);
param = param(:,mini);


%%% assign fitting parameters to the result struct 
tcfit.rmax = param(1);  tcfit.c50 = param(2);
tcfit.n = param(3);     tcfit.m = param(4);

%%% compute the explained varaince
r2 = 1 - ss / sum( (mn - mean(mn)).^2 );
if r2 <0; r2 = 0;end
if r2 >1; r2 = 1;end
tcfit.r2 = r2;

%%% raw tuning curve
tcfit.val = mn;
tcfit.co = co;


%%% fitted contrast response function
tcfit.x = 0:0.001:1;
tcfit.y = hyperratiofct(tcfit.x, tcfit.rmax, tcfit.c50, tcfit.n, tcfit.m);
tcfit.auc = sum(hyperratiofct( 0.1:0.1:1, tcfit.rmax, tcfit.c50, tcfit.n, tcfit.m));



if nargin==3
    %%% bootstrapping the fitting parameters
    parfor i = 1:1000
        bootidx = randi(length(mn), length(mn), 1); 
        boot(i).fitparam = fitCO(mn(bootidx), co_all(bootidx));
    end
    tcfit.boot = boot;
end


end



function ss = cf(x, y, param)
% cost function 
% x = contrast samples
% y = spike rates 
% param = parameters to fit the hyperbolic ratio function

rmax = param(1);    % maximal spike rate
c50 = param(2);     % contrast at half max response
n = param(3);       % exponent determining the steepness
m = param(4);       % tuning curve offset


y_pred = hyperratiofct( x, rmax, c50, n, m); % predicted data

ss = sum((y_pred-y).^2); %summed squared error, i.e. the general cost function 


% restrictions / boundaries
% the fit is restricted to positive paramaters and a saturation that has to
% be smaller than the maximal response
if c50 < 0 || n < 0 || m<0 || rmax <0 || rmax+m > max(y) 
    ss = Inf;
end
    
end


