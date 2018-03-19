function [ex, samplesize, c_all] = getPartialTrialSS(ex, ss)
%% limit ex.Trials to have the minimum repeats of 'ss'
% INPUT: ex
%             ss ... 8 maybe reasonable (used to be 4 for variability
%             anaylsis)
% OUTPUT: new ex
%                samplesize ... sample size [before after]
%
% written by Katsuhisa (18.07.17)
% +++++++++++++++++++++++++++++++++++

% initialization
samplesize = zeros(1,2);

% parse ex-file for the minimum repeat
stmtype = ex.exp.e1.type;
[u,c0] = uniquecount([ex.Trials.(stmtype)]);
samplesize(1) = sum(c0(1:end-1));
ko = zeros(1, length(ex.Trials));
for i = 1:length(ex.Trials)
    if ismember(ex.Trials(i).(stmtype), u(c0 < ss))
        ko(i) = 1;
    end
end
ex.Trials(ko==1) = [];
[~,c] = uniquecount([ex.Trials.(stmtype)]);
samplesize(2) = sum(c(1:end-1));

c_all = nan(1, length(u));
c_all(c0 >= ss) = c; 

