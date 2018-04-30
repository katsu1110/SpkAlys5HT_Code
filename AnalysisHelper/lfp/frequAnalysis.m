function ex = frequAnalysis( ex, varargin )
% filters and interpolates lfp traces and returns their the ex file also
% containing the filtered, interpolated and power spectra
%
% optional input variables are:
% fs : sampling frequency
% notchf : notch filtered frequency
% notchord : notch filter order
% bpf:   band pass filtered frequency
% bpord: band pass filter order
% nw    : time bandwidth product
%
% @CL 
%


%% define variables
stimdur = getStimDur(ex); % stimulus presentation duration

% notch filter variables
Fs = 1000;              % sampling frequency
% notchf = [49 51];     % notch filter frequency1
% notchf2 = [99 101];
% notchf3 = [149 151];
% notchf4 = [33 35];
% notchf5 = [66 68];
% notchord = 2;           % filter order

% bpf = [0.5 100];     % bandpass filter cutoff frequency
bpord = 2;            % filter order

% % power spectrum in Chalk et al. time halfbandwidth was 3 and k was 5
% nw = 2;         % time half bandwidth product
% 
% nfft = 1024; % number of frequency samples - I think

% time relative to stimulus onset to filter and interpolate LFP
t_off = -0.2;

% parameters for chronux toolbox
params = define_params;

% TODO: how to handle files with missing LFP entries (such as ma0014,
% 5.44pm)
if any(cellfun(@isempty, {ex.Trials.LFP}))
%     warning('TODO: in frequAnalysis line 44. How to proceed with empty LFP entries.');
    ex.Trials = ex.Trials(~cellfun(@isempty, {ex.Trials.LFP}));
end

% %% parse input
% k = 1;
% while k<=length(varargin)
%     switch varargin{k}
%         case 'fs'
%             Fs = varargin{k+1};
%         case 'notchf'
%             notchf = varargin{k+1};
%         case 'notchord'
%             notchord = varargin{k+1};
%         case 'bpf'
%             bpf = varargin{k+1};
%         case 'bpord'
%             bpord = varargin{k+1};
%         case 'nw'
%             nw = varargin{k+1};
%     end
%     k=k+2;
% end
% clearvars k;

%% generate filters
% define notch filter
% [b_notch,a_notch] = butter(notchord, notchf/(Fs/2), 'stop' );
% [b_notch2,a_notch2] = butter(bpord, notchf2/(Fs/2), 'stop');
% [b_notch3,a_notch3] = butter(bpord, notchf3/(Fs/2), 'stop');
d = designfilt('bandstopiir','FilterOrder', bpord, ...
               'HalfPowerFrequency1',49.5,'HalfPowerFrequency2',50.5, ...
               'DesignMethod','butter','SampleRate',Fs);
           
% % define bandpass filter
% [b_bpass,a_bpass] = butter(bpord, bpf/(Fs/2), 'bandpass');
% [b_bpass_notch,a_bpass_notch] = butter(bpord/2, notchf/(Fs/2), 'bandpass');

%% perform functions on each trial lfp
N = length(ex.Trials);
on = zeros(1, N); % overall noise
v = [];
for ind = 1:N
    % timing information
    t_frame = ex.Trials(ind).Start - ex.Trials(ind).TrialStart; % time of frame onsets
    t_lfp = ex.Trials(ind).LFP_ts - t_frame(1) ; % time rel:stimulus onset
    time = t_frame(1)+t_off : 1/Fs : t_frame(1)+stimdur;
    time = time - t_frame(1);

    % bandpass filter
%     filt1 = filtfilt(b_bpass, a_bpass, ex.Trials(ind).LFP);        
        filt1 = ex.Trials(ind).LFP;
        
    %% notch filter - filters line noise    
    % apply notch filter
%     filt2 = filtfilt(b_notch, a_notch, filt1);
%     filt2 = filtfilt(b_notch2, a_notch2, filt2);
%     filt2 = filtfilt(b_notch3, a_notch3, filt2);
    filt2 = filtfilt(d, filt1);
    
%     % remove 50Hz line noise with regression
%     filt2 = rmlinesc(filt2,params,0.05/N,0,50);

%     filt2 = ex.Trials(ind).LFP; % for debugging only
    % detrending
%     filt2 = locdetrend(filt2, Fs);

    % reduce the lfp signal to the period of stimulus presentation
    LFP_proc = interp1( t_lfp, filt2, time );
    
    % we are only interested in the stimulus induced fluctuations
%     LFP_proc = LFP_proc - nanmean(LFP_proc(time <= 0));
    v = [v; LFP_proc];
    
    % mutli taper function
    [ex.Trials(ind).POW, ex.Trials(ind).FREQ] = ...
        mtspectrumc(LFP_proc, params);
    on(ind) = mean(ex.Trials(ind).POW);
     
    % the preprocessed lfp and the corresponding time vector
    ex.Trials(ind).LFP_prepro = LFP_proc; 
    ex.Trials(ind).LFP_prepro_time = time;
end

ex.time = time; % time corresponding to saved lfp signal
ex.period = {[-0.2, 0], [0.05, 0.3], [0.35, time(end)]};

% remove trials with too much noise
ex.Trials(on > 3*std(on)) = [];
v(on > 3*std(on), :) = [];
 
% z-scoring
me = mean(v(:));
sd = std(v(:));
for ind = 1:length(ex.Trials)
    % z-scoring
    ex.Trials(ind).LFP_z = (ex.Trials(ind).LFP_prepro - me)/sd;
    
    % baseline, stimulus evoked, sustained
    for u = 1:3
        % LFP trace
        ex.Trials(ind).period(u).LFP_z_time = time(time>=ex.period{u}(1) & time <= ex.period{u}(2));
        ex.Trials(ind).period(u).LFP_z = ex.Trials(ind).LFP_z(time>=ex.period{u}(1) & time <= ex.period{u}(2));
         
        % multi taper function
        [ex.Trials(ind).period(u).POW, ex.Trials(ind).period(u).FREQ] = ...
            mtspectrumc(ex.Trials(ind).period(u).LFP_z, params);
    end
end


function avgstimdur = getStimDur(ex)
% returns the averaged and rounded stimulus presentation duration across
% trials


t_frame = cellfun(@(x, y) x-y, {ex.Trials.Start}, {ex.Trials.TrialStart}, ...
    'UniformOutput', 0); % time of frame onsets
stimdur = cellfun(@(x) x(end)+mean(diff(x)) - x(1), t_frame);


% also round the average stimulus duration to 2 digits precision
avgstimdur = round(mean(stimdur), 2);
