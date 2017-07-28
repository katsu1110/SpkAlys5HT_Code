function exinfo = bootstrap_exinfo( exinfo )
% evaluates the resampled tuning curves and returns the relative rate

if exinfo.isRC
    exinfo.pmodulation = 2;
    exinfo.resampled_mntc_base = 1;
    exinfo.resampled_mntc_drug = 1;
    exinfo.mntc_CI_base = 1; 
    exinfo.mntc_CI_drug = 1;
    return;
end

% get stimulus conditions to match indices
s1 = exinfo.ratepar;
s2 = exinfo.ratepar_drug;

% compute mean and resampled data
[rate_base, resampled_mntc_base] = ...
    getMnAndResDist(exinfo.ratemn, exinfo.rate_resmpl, exinfo.nrep, s1, s2);
[rate_drug, resampled_mntc_drug] = ...
    getMnAndResDist(exinfo.ratemn_drug, exinfo.rate_resmpl_drug, exinfo.nrep_drug, s2, s1);


exinfo.nonparam_ratio = mean(rate_drug)/mean(rate_base);

% return the confidence intervals of the averaged, resampled tuning curve
exinfo.mntc_CI_base = getCI(resampled_mntc_base); 
exinfo.mntc_CI_drug = getCI(resampled_mntc_drug); 

exinfo.resampled_mntc_base = resampled_mntc_base;
exinfo.resampled_mntc_drug = resampled_mntc_drug;


end


function [mnrates, mntc_resampled] = getMnAndResDist(mnrates, resmpls, nrep, i1, i2)
%returns the average response across stimulus conditions (i.e. the area
% under the tuning curve divided by the number of tested stimuli) 
% and the same metric for the bootstrapped responses.
%
%
% i1 and i2 are the stimuli in the baseline and drug condition. Note that
% it is important that we consider the same stimulus set when comparing the
% two experiments.


idx = ismember(i1, i2); % use only those stimuli that were used in both experiments

mnrates = mnrates(idx); % average spike rates (raw tc)
resmpls = resmpls(idx); % re-samples, cell array
% nrep = nrep(idx); %only needed for weighted average


% each column in each cell contains a set of resampled spiking responses
% compute the mean for each of these sets to obtain the distribution for
% each data point in the tuning curve
mntc_resampled = nan(length(resmpls), 1000);
for i = 1:length(resmpls)
    mntc_resampled(i,:) = mean(resmpls{i},1);
end

% now we can average across the resampled and averaged responses
try
%     if size(nrep, 1)>1
%         nrep = nrep';
%     end
%     res_mnmn = nrep*res_mn / sum(nrep); % weighted average
    mntc_resampled = mean(mntc_resampled); 
catch
    disp('');
end


end


function ci = getCI(y)
% 25, 50 and 75% intervals
% allows for Tukey's test of outliers

ci(1) = prctile(y, 5);
ci(2) = prctile(y, 25);
ci(3) = prctile(y, 50);
ci(4) = prctile(y, 75);
ci(5) = prctile(y, 95);


% % k times the quartiles as range for outliers
% k = 0;
% ci(1) = ci(1)-k*(ci(3)-ci(1));
% ci(2) = ci(2)-k*(ci(3)-ci(2));
% ci(4) = ci(4)+k*(ci(4)-ci(3));
% ci(5) = ci(5)+k*(ci(5)-ci(3));

end


