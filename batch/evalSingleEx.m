function [ex, argout] = evalSingleEx( exinfo, fname, varargin)
% [ex, argout] = evalSingleEx( exinfo, fname )
% 
% This function is only called by runExinfoAnalysis. Dependign on the
% experiment (drifting or flashed grating), it calls the function
% performing either  reverse correlation analysis or classic spiking
% analysis. Arguments are passed on. See their documentation for more
% detail.
%
% 
% @CL


if exinfo.isRC
    [ex, argout] = evalSingleRC(exinfo, fname, varargin{:});
else
    [ex, argout] = evalSingleDG(exinfo, fname, varargin{:});
end


end



