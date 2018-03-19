function countSampleSize(dat)
%% simply count the sample size for the noise correlation analysis
%
% written by Katsuhisa (19.07.17)
% +++++++++++++++++++++++++++++++++++++++++++++++

c = 0;
num = length(dat.expInfo);
for i = 1:num
    c = c + sum(dat.expInfo(3).ff.classic_2ndhalf.stimrep(1:end-1));
end
disp(['total sample size: ' num2str(c)])
disp(['the number of experiment: ' num2str(num)])

% For noise correlation analysis
% OR: c = 4750, N = 95
% SF: c = 1188, N = 36
% CO: c = 3375, N = 125
% SZ: c = 1824, N = 48