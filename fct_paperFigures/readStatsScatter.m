function [n5HT, nNaCl, pwil, ptt, r_p, r_s] = readStatsScatter(s)
% reads out the statistical values from the scatter plot title
% s is the titles's string.
%
% @CL 15.12.2016

% #5HT
[~,i_strt] = regexp(s{1}, '(n='); i_end = regexp(s{1}, ')');
n5HT = str2double(s{1}(i_strt+1:i_end-1));

% #NaCl
[~,i_strt] = regexp(s{2}, '(n='); i_end = regexp(s{2}, ')');
nNaCl = str2double(s{2}(i_strt+1:i_end-1));

% p value paired t-test
[~, i_strt] = regexp(s{3}, 'paired ttest p='); 
ptt = str2double(s{3}(i_strt+1:i_strt+5));

% p value paired, signrank test
[~, i_strt] = regexpi(s{3}, 'signrank test p='); 
pwil = str2double(s{3}(i_strt+1:i_strt+5));



% correlation values (5HT only)
[~, i_strt] = regexpi(s{1}, 'r_{p}='); 
r_p = getR(s{1}, i_strt);

% correlation values (5HT only)
[~, i_strt] = regexpi(s{1}, 'r_{s}=');
r_s = getR(s{1}, i_strt);

end



function r = getR(s, i_strt)

if strcmp(s(i_strt+1), '-')
    r(1) = str2double(s(i_strt+1:i_strt+5));
    r(2) = str2double(s(i_strt+9:i_strt+14));
else
    r(1) = str2double(s(i_strt+1:i_strt+4));
    r(2) = str2double(s(i_strt+8:i_strt+13));
end

end

