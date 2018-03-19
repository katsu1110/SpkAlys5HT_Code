function [ex_sps, ex_lps, ex] = pupilSplit(ex)
%% split ex.Trials by median split of the mean pupil size
%
% OUTPUT: ex_sps ... small pupil size
%                ex_lps ... large pupil size
% 
% written by Katsuhisa (28.09.17)
% +++++++++++++++++++++++++++++++++++++

% label the trials
[label_seq] = label4StmPerTr(ex);

% old time or new time of ex-file structure
if isfield(ex.Trials(1),'TrialStart_remappedGetSecs')
      time = 'N';   % new
elseif ~isfield(ex.Trials(1),'TrialStart_remappedGetSecs')...
        && isfield(ex.Trials(1),'TrialStartDatapixx')...
        && isfield(ex.Trials(1).times, 'lastEyeReadingDatapixx')
      time = 'O';   % old
end

% raw pupil size
len_tr = length(ex.Trials);
ps = [];
for i = 1:len_tr
    % n-th stimulus
    ex.Trials(i).n_stm = label_seq(i);

    % timing of start and end of stimulus presentation
    if strcmp(time, 'N')            
        t = ex.Trials(i).Eye.t(1:ex.Trials(i).Eye.n) - ex.Trials(i).TrialStartDatapixx;
        st = ex.Trials(i).Start - ex.Trials(i).TrialStart_remappedGetSecs;           
    elseif time == 'O'
        delta = ex.Trials(tr(i)).Eye.t(ex.Trials(tr(i)).Eye.n) - ex.Trials(tr(i)).TrialStartDatapixx - ex.Trials(tr(i)).times.lastEyeReading;
        t = ex.Trials(tr(i)).Eye.t(1:ex.Trials(tr(i)).Eye.n)-ex.Trials(tr(i)).TrialStartDatapixx-delta;
        st = ex.Trials(tr(i)).Start - ex.Trials(tr(i)).TrialStart;
    end

    % get the timing of start and end of stimulus
    [~,stpos] = min(abs(t-st(1)));
    [~,enpos] = min(abs(t-st(end)));         

    temp = nanmedian([ex.Trials(i).Eye.v(3,:); ex.Trials(i).Eye.v(6,:)],1);
    temp = temp(~isnan(temp) & ~isinf(temp));    
    ex.Trials(i).pupil_raw = temp(stpos:enpos);
    ps = [ps ex.Trials(i).pupil_raw];
end

me = nanmean(ps);
sd = nanstd(ps);

% z-scoring and extract pupil value
for i = 1:len_tr              
        % z-scoring
        ex.Trials(i).pupil_raw = (ex.Trials(i).pupil_raw - me)/sd;    

        % last 1/4
        l = length(ex.Trials(i).pupil_raw);
        ex.Trials(i).pupil_val = nanmean(ex.Trials(i).pupil_raw(end-round(l/4)+1:end));

%             % mean of ps
%             ex.Trials(i).pupil_val = nanmean(ex.Trials(i).pupil_raw);
end

% median split based on the pupil_val on each stimulus type
[ stimparam, vals] = getStimParam( ex );
sps_tr = [];
lps_tr = [];
for s = 1:length(vals)
    idx = find([ex.Trials.(stimparam)] == vals(s));
    med = median([ex.Trials(idx).pupil_val]);
    sps_tr = [sps_tr, idx([ex.Trials(idx).pupil_val] < med)];
    lps_tr = [lps_tr, idx([ex.Trials(idx).pupil_val] > med)];
end
ex_sps = ex;
ex_sps.Trials = ex.Trials(sps_tr);
ex_lps = ex;
ex_lps.Trials = ex.Trials(lps_tr);
    
