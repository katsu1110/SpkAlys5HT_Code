function f1f0 = getPhaseSelectivity( ex, varargin )
% returns the well established f1 over f0 in an experiment to categroize
% simple and complex cells
% (f1/f0>1 simple cell, f1/f0<1 complex cells)
%
% @Corinna Lorenz, 22.08.2016


p_flag = false;     stim = 'tf';
k=1;
while k<=length(varargin)
    switch varargin{k}
        case 'plot'
            p_flag = varargin{k+1};
        case 'stim'
            stim = varargin{k+1};
            
    end
    k=k+1;
end



%%
% sdfs kernel
kern = ones(40,1)/40;

% spontaneoue activtiy from blank trials
blanktrials = ex.Trials([ex.Trials.Reward]==1 & [ex.Trials.(stim)]>1000);
if ~isempty(blanktrials)
    psth_blank = getPSTH(blanktrials); spontresp = mean(psth_blank);
else
    spontresp = 0;
end

% correct trials only
ex.Trials = ex.Trials([ex.Trials.Reward]==1 & [ex.Trials.(stim)]<1000);

allstim = unique([ex.Trials([ex.Trials.(stim)]<1000).(stim)]);
if strcmp(stim, 'tf')
    tf =  allstim;
else
    tf = repmat(ex.stim.vals.tf, size(allstim));
end

% f1 from fourier transformation and f0 as mean - spont response
for k = 1:length(allstim)
    
    %1.first derive the psth for all trials with particular stimuli
    trials = ex.Trials( [ex.Trials.(stim)] == allstim(k));
    
    [psth, t] = getPSTH(trials);
    psth = psth - spontresp;    psth(psth<0) = 0;
    
    try
        sdfs(k,1:length(psth)) = filter(kern,1, psth);

        % neglect the first 50ms corresponding to the general latency
        x = t(50:end);
        y = sdfs(k,50:length(psth)) ;
        [a(k),f0(k)] = f1f0Sine (y, x, tf(k));
                
        % alternatively compute the fft
        [pow, freq] = pmtm(y, 3, length(y)*10, 1000);
        
        freq_new = 0:0.01:tf(k)+0.01;
        pow_new = interp1(freq, pow, freq_new, 'spline');
        f0_(k) = mean(y);
        f1_(k) = pow_new( find(freq_new >= tf(k), 1, 'first'));
        
    catch
        disp('')
    end
    
end

t = 0.001:0.001:length(sdfs)/1000;
f1 = a; % f1 is the amplitude

% use the highest mean response to compute
[~, i] = max(f1);
f1f0 = f1(i)/f0(i);


if p_flag
    %%% plot results
    figure('Position', [   836   247   545   701]);
    
    offset = max(max(sdfs));
    for k = 1:length(allstim)
        subplot(length(allstim), 2, (k*2)-1)
        plot(t, sin(t *pi*2*tf(k) ) * offset/2+ offset/2, 'Color', [0.8 0.8 0.8]); hold on;
        plot(t, sdfs(k,:), 'b');
        title(num2str(allstim(k)));
        ylabel('spk/s');
        xlim([t(1) t(end)])
    end
    
    
   
    subplot(2, 2, 2);
    plot(allstim, f1./f0, 'x-');
    xlabel tf; ylabel f1/f0; xlim([allstim(1) allstim(end)]);
    title(sprintf('stim @maxresp: %1.1f, f1f0: %1.3f', allstim(i), f1f0));
    
    
    subplot(2,2,4)
    plot(allstim, f1, 'b'); hold on;
    plot(allstim, f0, 'b--'); legend('f1', 'f0');
    
    xlabel stim; xlim([allstim(1) allstim(end)]);
    
end
end


