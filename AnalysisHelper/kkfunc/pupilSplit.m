function [ex_sps, ex_lps, ex] = pupilSplit(ex)
%% preprocess pupil size and split ex.Trials by median split of the mean pupil size
%
% OUTPUT: ex_sps ... small pupil size
%                ex_lps ... large pupil size
% 
% written by Katsuhisa (28.09.17)
% +++++++++++++++++++++++++++++++++++++

% label the trials
[label_seq] = label4StmPerTr(ex);
isRC = 1;
if max(label_seq) > 1
    isRC = 0;
end

% completed trials
completed = find([ex.Trials.Reward]>0);

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
psme = [];
for i = 1:len_tr
    % n-th stimulus
    ex.Trials(i).n_stm = label_seq(i);
    ex.Trials(i).pupil_raw = nan;
    if ismember(i, completed)
        % timing of start and end of stimulus presentation
        if strcmp(time, 'N')            
            t = ex.Trials(i).Eye.t(1:ex.Trials(i).Eye.n) - ex.Trials(i).TrialStartDatapixx;
            st = ex.Trials(i).Start - ex.Trials(i).TrialStart_remappedGetSecs;           
        elseif time == 'O'
            delta = ex.Trials(i).Eye.t(ex.Trials(i).Eye.n) - ex.Trials(i).TrialStartDatapixx - ex.Trials(i).times.lastEyeReading;
            t = ex.Trials(i).Eye.t(1:ex.Trials(i).Eye.n)-ex.Trials(i).TrialStartDatapixx-delta;
            st = ex.Trials(i).Start - ex.Trials(i).TrialStart;
        end

        % get the timing of start and end of stimulus
        [~,stpos] = min(abs(t-st(1)));
        [~,enpos] = min(abs(t-st(end)));         

        % pupil data
        temp = nanmedian([ex.Trials(i).Eye.v(3,:); ex.Trials(i).Eye.v(6,:)],1);
        temp = temp(~isnan(temp) & ~isinf(temp));    
        ex.Trials(i).pupil_raw = temp(stpos:enpos);    
        psme = [psme, nanmean(ex.Trials(i).pupil_raw)];
    end
end

% filter the pupil size
ps = [];
[mseq] = HiPaFi(psme);      % prepare high-pass filter (my function)
c = 1;
for i = 1:len_tr
    ex.Trials(i).pupil_filt = nan;
    if ismember(i, completed)
        % band-pass filtering
        ex.Trials(i).pupil_filt = ex.Trials(i).pupil_raw - mseq(c);                    % high-pass filtering
        ex.Trials(i).pupil_filt = LoPaFi(ex.Trials(i).pupil_filt, 1);      % low-pass filtering (my function) 
        ps = [ps ex.Trials(i).pupil_filt];
        c = c + 1;
    end
end

% z-scoring and extract pupil value
me = nanmean(ps);
sd = nanstd(ps);
for i = 1:len_tr             
    ex.Trials(i).pupil_z = nan;
    if ismember(i, completed)
        % z-scoring
        ex.Trials(i).pupil_z = (ex.Trials(i).pupil_filt - me)/sd;    
    end
end
% pupil metric
for i = 1:len_tr  
    ex.Trials(i).pupil_val = nan;
    if ismember(i, completed)
        if isRC==1 % last 1/8 (250ms) --- RC
            l = length(ex.Trials(i).pupil_z);
            ex.Trials(i).pupil_val = nanmean(ex.Trials(i).pupil_z(end-round(l/8)+1:end));
        else % 4 stimuli per trial
            try
                if label_seq(i + 4 - label_seq(i))==4
                    ex.Trials(i).pupil_val = nanmean(ex.Trials(i+4-label_seq(i)).pupil_z);
                end
            catch
                  continue;
            end
        end
    end
end
% remove trials with pupil_val = nan
ex.Trials = ex.Trials(~isnan([ex.Trials.pupil_val]));
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
    
