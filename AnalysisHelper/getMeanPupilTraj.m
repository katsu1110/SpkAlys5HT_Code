function v3 = getMeanPupilTraj( ex )
% retracts the pupil size over time 

 
% interpolate all v3 to similar structure
v3 = cellfun( @getV3, ...
    {ex.Trials.Eye}, {ex.Trials.TrialStartDatapixx}, ...
    {ex.Trials.TrialEndDatapixx}, ...
    'UniformOutput', false);

mean(v3);

end



%-------------------------------------------------------------------------

function rsz_multiple = getRSZ(rsz_single, t)
rsz_multiple = zeros(size(t));
rsz_multiple(find(t>=0, 1, 'first')) = .1;
rsz_multiple(end) = rsz_single;
end

function v3 = getV3(eye, trialstart, trialend)
 
dur = trialend - trialstart;
t = eye.t(1:eye.n) - trialstart;
v3 = eye.v(3,1:eye.n);

i_pre = t<0;
i_post = t>trialend;
i_dur = ~i_pre & ~i_post;


t_pre_new = -0.02:0.001:0;
v3_pre = interp1(t(i_pre), v3(i_pre), t_pre_new);


% filter all v3
fc = 20; %cutoff frequency in Hz
[b,a] = butter(3, fc*2/500); %lowpass filter

dat.v3 = filtfilt(b,a, v3);



end

