function [para] = LFPbyStm(ex)
% compute LFP parameters in the specified ex-file, split by types of
% stimuli
% - ex file must have LFP data, using 'loadCluster.m'
% - ex file must be filtered by 'filterLFP.m'

% stimulus type ==============
[ stimparam, vals] = getStimParam( ex );
% blank_i =  vals > 1;
lenv = length(vals);

% initialization ==================
comp = find(abs([ex.Trials.Reward])>0, 1, 'first');
fidx = ex.Trials(comp).FREQ>= 0 & ex.Trials(comp).FREQ<=100;
freq = ex.Trials(comp).FREQ(fidx);
fs = 1000;

% structure
para.stm.param = stimparam;
para.stm.vals = vals;
para.ts = ex.time;
para.f = freq';
para.window = ex.period;
para.ts = ex.time;
para.ntr = length(ex.Trials);
para.bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
para.range = {[0.4,3.9], [4, 7], [8, 13], [14, 29], [30, 80]};
lenb = length(para.bands);
para.params = define_params;
para.wnd = 0.1; % was 0.3
[~, zerot] = min(abs(-para.wnd:1/1000:para.wnd));
% ampwindow = wndrange > -wnd &  wndrange  < wnd;
for i = 1:lenv    
    % trials with the unique stimulus
    trials = ex.Trials([ex.Trials.(stimparam)] == vals(i));
    lentr = length(trials);
    
    % firing rate per trial
    [~, spkc]  = getSpks(trials, [0 0]);
    para.stm.fr(i) = sum(spkc)/(lentr/ex.fix.duration);

    % LFP trace averaged trials across the same stimulus
    lfpfull = vertcat(trials.LFP_z);
    lfpfull = lfpfull(mean(isnan(lfpfull), 2)==0, :);
    para.lfp_stm.mean(i,:) = mean(lfpfull, 1);
    para.lfp_stm.sem(i,:) = std(lfpfull, [], 1)/sqrt(lentr);

    % placeholder
    for u = 1:3
        ph.period(u).stlfp.mat = [];
        ph.period(u).lfpz = [];
        for b = 1:lenb
%             para.period(u).(['lfp_' para.bands{b} '_tc']) = [];
            ph.period(u).(['lfp_' para.bands{b} '_pow']) = [];
        end
        ph.period(u).pow = [];
    end
    %      % CL's correction
%        para.stlfp.avg_stlfp(i,:) = para.stlfp.avg_stlfp(i,:) ...
%            - mean(para.stlfp.avg_stlfp(i, wndrange < -0.06));
    for n = 1:lentr
        for u = 1:3
            % lfp trace
            ph.period(u).lfpz = [ph.period(u).lfpz; trials(n).period(u).LFP_z]; 
            % spike-triggered average LFP
            spk = getSpks(trials(n), [para.window{u}(1), 2-para.window{u}(2)]);
            stlfp = getSTA(trials(n).period(u).LFP_z, trials(n).period(u).LFP_z_time, ...
                    spk{1}, para.wnd, fs);                
            ph.period(u).stlfp.mat = [ph.period(u).stlfp.mat; stlfp];
            % frequency band
            for b = 1:lenb
%                 para.period(u).(['lfp_' para.bands{b} '_tc']) = [para.period(u).(['lfp_' para.bands{b} '_tc']); ...
%                     trials(n).period(u).(['lfp_' para.bands{b} '_tc'])];
                ph.period(u).(['lfp_' para.bands{b} '_pow']) = [ph.period(u).(['lfp_' para.bands{b} '_pow']); ...
                    trials(n).period(u).(['lfp_' para.bands{b} '_pow'])];
            end
            % averaged power
            ph.period(u).pow = [ph.period(u).pow; trials(n).period(u).POW'];
        end
    end
    for u = 1:3
        ph.period(u).lfpz = ph.period(u).lfpz(any(isnan(ph.period(u).lfpz), 2)==0, :);
    end
    
     % baseline, stimulus evoked, sustained
    for u = 1:3
        % spike-triggered LFP
        para.period(u).stlfp.nspk(i) = size(ph.period(u).stlfp.mat, 1);
        para.period(u).stlfp.avg_stlfp(i,:) = sum(ph.period(u).stlfp.mat, 1)/para.period(u).stlfp.nspk(i);
        para.period(u).stlfp.sem_stlfp(i,:) = nanstd(ph.period(u).stlfp.mat, [], 1)...
            /sqrt(para.period(u).stlfp.nspk(i));
        
        % STA baseline correction
        para.period(u).stlfp.avg_stlfp(i,:) = para.period(u).stlfp.avg_stlfp(i,:)...
            - mean(para.period(u).stlfp.avg_stlfp(i,1:floor(para.wnd*fs/5)), 2);
        
        para.period(u).stlfp.zeroamp(i) = para.period(u).stlfp.avg_stlfp(i, zerot);
        [para.period(u).stlfp.POW(:,i), para.period(u).stlfp.FREQ(:,i)] =...
            mtspectrumc(para.period(u).stlfp.avg_stlfp(i,:), para.params);
        f = para.period(u).stlfp.FREQ(:,i);
        for b = 1:lenb
            para.period(u).stlfp.(['pow_' para.bands{b}]) = ...
                nanmean(para.period(u).stlfp.POW(f >= para.range{b}(1) & f <= para.range{b}(2), i));
        end
        
        % spectrogram
        [para.period(u).spectrogram.S{i}, para.period(u).spectrogram.t{i}, para.period(u).spectrogram.f{i}] = ...
            mtspecgramc(ph.period(u).lfpz', [0.5, 0.05], para.params);
%         [0.5 0.05]
%         t{i} = t{i} + ex.Trials(end).LFP_prepro_time(1);
        
        % spike-LFP coherency
        spk = getSpks(trials, [para.window{u}(1), 2 - para.window{u}(2)]); 
%         spk = spk(mean(isnan(para.period(u).lfpz),2)==0);
        [para.period(u).coherence.C{i}, para.period(u).coherence.phi{i}, ...
            para.period(u).coherence.S12{i}, para.period(u).coherence.S1{i}, ...
            para.period(u).coherence.S2{i}, para.period(u).coherence.f{i}]= ...
            coherencycpt(ph.period(u).lfpz', cell2struct(spk, 'spk', 1)', para.params);

        % frequency band
        for b = 1:lenb
%             para.period(u).lfp_stm_wave(b).mean(i,:) = nanmean(para.period(u).(['lfp_' para.bands{b} '_tc']), 1);
            para.period(u).lfp_stm_wave(b).pow(i) = nanmean([ph.period(u).(['lfp_' para.bands{b} '_pow'])], 1);
        end
        
        % averaged power across the same stimulus
        para.period(u).pow_avg(i,:) = nanmean(ph.period(u).pow, 1);    
        para.period(u).freq(i,:) = trials(1).period(u).FREQ;
    end
end