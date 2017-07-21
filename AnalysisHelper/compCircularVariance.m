function v = compCircularVariance( r, theta )
% v = compCircularVariance( R, theta )
%
% compute the circular variance according to Ringach et al. (2002)
% with 
%   R = sum(r e^(i2theta)) / sum(r)
%   V = 1-abs(R)
%
% where theta are the stimulus orientations and r the corresponding
% responses from a recorded unit.
%
%
%
% see here for the original paper online
% http://www.jneurosci.org/content/22/13/5639.long
%
%
% @CorinnaLorenz, 11/07/2017



k = exp(1i*2*theta);
R = sum(r.*k)/sum(r);

v = 1-abs(R);

end

