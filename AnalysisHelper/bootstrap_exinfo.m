function exinfo = bootstrap_exinfo( exinfo, conc_flag )
% evaluates the resampled tuning curves and returns the relative rate

if exinfo.isRC
    return;
end

% get stimulus conditions to match indices
s1 = exinfo.ratepar;
s2 = exinfo.ratepar_drug;

% compute mean and resampled data
[rate_base, res_tc_base] = ...
    getMnAndResDist(exinfo.ratemn, exinfo.rate_resmpl, exinfo.nrep, s1, s2);
[rate_drug, res_tc_drug] = ...
    getMnAndResDist(exinfo.ratemn_drug, exinfo.rate_resmpl_drug, exinfo.nrep_drug, s2, s1);


exinfo.nonparam_ratio = mean(rate_drug)/mean(rate_base);

% return the confidence intervals of the averaged, resampled tuning curve
exinfo.mntc_CI_base = getCI(res_tc_base); 
exinfo.mntc_CI_drug = getCI(res_tc_drug); 

if nargin == 2 && conc_flag
    exinfo.rate_resmpl = res_tc_base;
    exinfo.rate_resmpl_drug = res_tc_drug;
end
end


function [ratemn, res_mnmn] = getMnAndResDist(ratemn, resmpls, nrep, i1, i2)
%returns data with similar rate 

idx = ismember(i1, i2);
ratemn = ratemn(idx); % mean spike rate (across raw tc)
resmpls = resmpls(idx); % resamples, cell array
nrep = nrep(idx);

% compute the average tuning curve from the resampled stimulus conditions
res_mn = nan(length(resmpls), 1000);
for i = 1:length(resmpls)
    res_mn(i,:) = mean(resmpls{i},1);
end

try
    if size(nrep, 1)>1
        nrep = nrep';
    end
    res_mnmn = nrep*res_mn / sum(nrep); % weighted average

catch
    disp('');
end
% 
% resmpls = resample2(ratemn);
% res_mnmn = mean(res_mn, 1);

end


function ci = getCI(y)
% 25, 50 and 75% intervals
% allows for Tukey's test of outliers
ci(1) = prctile(y, 5);
ci(2) = prctile(y, 25);
ci(3) = prctile(y, 50);
ci(4) = prctile(y, 75);
ci(5) = prctile(y, 95);

% k times the quartiles as range for outliers
k = 0;

ci(1) = ci(1)-k*(ci(3)-ci(1));
ci(2) = ci(2)-k*(ci(3)-ci(2));
ci(4) = ci(4)+k*(ci(4)-ci(3));
ci(5) = ci(5)+k*(ci(5)-ci(3));

end


function res = resample2(A)
% resample from A

n = length(A); 
nsmpl = 1000; % number 

res = ones(n, nsmpl); %initial variable to prevent overhead
idx = randi(n, 1000, n); % combination of indices corresponding to A's data

for i = 1:nsmpl
    res(:, i) = A(idx(i, :)); %assigning the randomized resamples
end

end


