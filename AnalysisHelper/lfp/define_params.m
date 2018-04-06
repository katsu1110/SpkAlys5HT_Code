function params = define_params
% define parameter sets for chronux toolbox

params.err = 0;
params.Fs = 1000;
params.fpass = [0 100];
params.tapers = [3,5]; % was [2,3] but Chalk et al.,2010 used [3,5]
params.pad = 0;
params.trialave = 1;