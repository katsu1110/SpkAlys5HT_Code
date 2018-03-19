function [ex, spkrate] = znormex(ex, exinfo, pupil, varargin)
%[ex, spkrate, spkcount] = znormex(ex, exinfo)
% 
% computes the stimulus elicited spike count and resulting spike rate and
% saves the z-normed activity in ex.Trials  for each stimulus.
% 
% If 'pupil' is 0, not median split by the pupil value was performed. 
%
%
% additionally, cell arrays with information for the raster plot are saved
% in ex (ex.raster, ex.rastercum, ex.raster, ex.trials_n). The cells are
% sorted corresponding to spkrate.x with x being the stimulus feature (or,
% co, sz, sf,...)
% 
% @CL
%
% pupil analysis added by Katsuhisa (06.06.17)
% ++++++++++++++++++++++++++++++++++++++++++++++++++

% preprocess pupil
if nargin < 3
        pupil = 0;
end

if pupil==1
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
    for i = 1:length(ex.Trials)                
            % z-scoring
            ex.Trials(i).pupil_raw = (ex.Trials(i).pupil_raw - me)/sd;    

            % last 1/4
            l = length(ex.Trials(i).pupil_raw);
            ex.Trials(i).pupil_val = nanmean(ex.Trials(i).pupil_raw(end-round(l/4)+1:end));

%             % mean of ps
%             ex.Trials(i).pupil_val = nanmean(ex.Trials(i).pupil_raw);
    end
    
    ps_s_ind = [];
    ps_l_ind = [];
end

        
% stimulus
param1 = exinfo.param1;
parvls = unique( [ ex.Trials.(param1) ] );

j = 1;
Fs = 1000; % recording frequency

if exinfo.isadapt
    time = 0:1/Fs:5;
else
    time = 0.001:1/Fs:0.45;
end

for par = parvls
    
    % z-scored spikes for this stimulus
    ind = par == [ ex.Trials.(param1) ];
    spkrate(j).nrep = sum(ind);
    spkrate(j).(param1) = par;

    
    % spike rates statistics
    spkrate(j).mn  = mean([ex.Trials(ind).spkRate]);
    spkrate(j).var = var([ex.Trials(ind).spkRate]);
    spkrate(j).sd = std([ex.Trials(ind).spkRate]);
    spkrate(j).sem = spkrate(j).sd / sqrt( spkrate(j).nrep );
    spkrate(j).raw = [ex.Trials(ind).spkRate];
    spkrate(j).resamples = resample2([ex.Trials(ind).spkRate]);
  
    % spike count statistics
    spkrate(j).cmn  = mean([ex.Trials(ind).spkCount]);
    spkrate(j).cvar = var([ex.Trials(ind).spkCount]);
    
    index = find(ind);
    
    if pupil==1 
            % median split by pupil for spike rate
            med = median([ex.Trials(ind).pupil_val]);

            % small pupil trials
            sind = index([ex.Trials(ind).pupil_val] < med);        
%             spkrate(j).smallps = sind;
            ps_s_ind = [ps_s_ind sind];

            % large pupil trials
            sind = index([ex.Trials(ind).pupil_val] > med);          
%             spkrate(j).largeps = sind;
            ps_l_ind = [ps_l_ind sind];
    end
       
    % z normed spike counts
    z = zscore( [ ex.Trials(ind).spkCount ] );
    z = num2cell(z);
    [ex.Trials(ind).zspkcount] = deal(z{:});
    
    % raster - contains 0 for times without spike and 1 for times of spikes
    % is reduced for times between -0.5s and 0.5s around stimulus onset
    idxct = find(ind);
    ct = ex.Trials(idxct);    % current trials
    
    ex.trials_n(j) = sum(ind);
    ex.raster{j,1} = zeros(length(ct), length(time));
    ex.rastercum{j,1} = nan(length(ct), length(time));
        
    % convert the spike times into a matrix raster with bits
    for k = 1:length(idxct)
        
        t_strt = ct(k).Start - ct(k).TrialStart;
        if exinfo.isadapt
            t_strt = t_strt(t_strt > ex.stim.vals.adaptationDur);
        end
        
        idx = round((ct(k).Spikes(...
            ct(k).Spikes>=t_strt(1) & ct(k).Spikes<=t_strt(end))-t_strt(1)) ...
            *Fs);
        
        idx(idx==0) = 1; % avoid bad indexing
        
        ex.raster{j}(k, idx) = 1;
        ex.rastercum{j}(k, idx) = k;
        
    end
    
    j = j+1;
    
end

if pupil==1
        % split ex-file based on pupil size
        ex.pupil_split(1).idx = ps_s_ind;
        ex.pupil_split(2).idx = ps_l_ind;
end

ex.raster_pval = parvls;

end



function res = resample2(A)
% resample from A

n = length(A); 
nsmpl = 1000; % number 

res = ones(n, nsmpl); %initial variable to prevent overhead
try
        idx = randi(n, nsmpl, n); % combination of indices corresponding to A's data
        for i = 1:nsmpl
            res(:, i) = A(idx(i, :)); %assigning the randomized resamples
        end
catch
        res = [];
end

end








