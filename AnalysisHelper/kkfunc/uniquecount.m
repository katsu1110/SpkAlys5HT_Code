function [u,c] = uniquecount(v)
%% count how many times unique values appear in the given vector v
%
%
% written by Katsuhisa (14.06.17)
% +++++++++++++++++++++++++++++++++++++++++++++++++

u = unique(v);
c = u;
for k = 1:length(u)
    c(k) = sum(v==u(k));
end