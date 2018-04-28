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
                    
% structure
para.stm.param = stimparam;
para.stm.vals = vals;
para.ts = ex.time;
para.f = freq';
para.period = ex.period;
para.ts = ex.time;
para.ntr = length(ex.Trials);
para.bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
para.range = {[0,4], [4, 7], [8, 13], [14, 29], [30, 80]};
lenb = length(para.bands);
para.params = define_params;

wnd = 0.1; % was 0.3
[~, zerot] = min(abs(-wnd:1/1000:wnd));
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
    para.lfp_stm.se(i,:) = std(lfpfull, [], 1)/sqrt(lentr);

     % baseline, stimulus evoked, sustained
    for u = 1:3
        peri = trials.period(u);
        
        % spectrogram
        lfp = vertcat(peri.LFP_z);
        lfp = lfp(mean(isnan(lfp), 2)==0, :);
        [para.period(u).spectrogram.S{i}, para.period(u).spectrogram.t{i}, para.period(u).spectrogram.f{i}] = ...
            mtspecgramc(lfp, [0.5 0.05], para.params);
%         t{i} = t{i} + ex.Trials(end).LFP_prepro_time(1);
        
        % spike-LFP coherency
        spk = getSpks(trials, [ex.period(1), 2 - ex.period(2)]); 
        spk = spk(mean(isnan(lfp),2)==0);
        [para.period(u).coherence.C{i}, para.period(u).coherence.phi{i}, ...
            para.period(u).coherence.S12{i}, para.period(u).coherence.S1{i}, ...
            para.period(u).coherence.S2{i}, para.period(u).coherence.fg{i}]= ...
            coherencycpt(lfp, cell2struct(spk, 'spk', 1), para.params);

        % frequency band
        for b = 1:lenb
            para.period(u).lfp_stm_wave(b).mean(i,:) = nanmean(vertcat(peri.(['lfp_' para.bands{b} '_tc'])), 1);
            para.period(u).lfp_stm_wave(b).pow(i) = nanmean([peri.(['lfp_' para.bands{b} '_pow'])]);
        end
        
        % averaged power across the same stimulus
        para.period(u).pow_avg(i,:) = nanmean(horzcat(peri.POW),2)';

        % spike-triggered LFP
        para.period(u).stlfp.mat = [];
        %      % CL's correction
    %        para.stlfp.avg_stlfp(i,:) = para.stlfp.avg_stlfp(i,:) ...
    %            - mean(para.stlfp.avg_stlfp(i, wndrange < -0.06));
        
    end
    
    % spike-triggered LFP
    for n = 1:lentr
        for u = 1:3
            stlfp = getSTA(trials(n).period(u).LFP_z, trials(n).period(u).LFP_z_time, ...
                    trials(n).Spikes(time>=para.period{u}(1) & time <= para.period{u}(2)), wnd);                
            para.period(u).stlfp.mat = [para.period(u).stlfp.mat; stlfp];
        end
    end
    for u = 1:3
        para.period(u).stlfp.nspk = size(para.period(u).stlfp.mat, 1);
        para.period(u).avg_stlfp = sum(para.period(u).stlfp.mat, 1)/para.period(u).stlfp.nspk;
        para.period(u).sem_stlfp = nanstd(para.period(u).stlfp.mat, 1)...
            /sqrt(para.period(u).stlfp.nspk);
        para.period(u).stlfp.zeroamp = para.period(u).stlfp.avg_stlfp(zerot);
        [para.period(u).stlfp.POW, para.period(u).stlfp.FREQ] =...
            mtspectrumc(para.period(u).avg_stlfp, para.params);
        f = para.period(u).stlfp.FREQ;
        for b = 1:lenb
            para.period(u).stlfp.(['pow_' para.bands{b}]) = ...
                nanmean(para.period(u).stlfp.POW(f >= para.range{b}(1) & f <= para.range{b}(2)));
        end
    end
end