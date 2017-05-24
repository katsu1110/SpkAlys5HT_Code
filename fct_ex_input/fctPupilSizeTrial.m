function [dat, oTrials] = fctPupilSizeTrial(ex, fname, p_flag)
% extracts pupil size and corresponding position in 2 sec trial period
% assigns the data into the trial struct


if nargin < 3
    p_flag = 1;    
end

if p_flag
    h = figure;
end


oTrials = ex.Trials;

fc = 20; %cutoff frequency in Hz
[b,a] = butter(3, fc*2/500);


col = lines(4);
j = 1;

% window length
if isempty(strfind(fname, 'RC'))
    len = 224;
else 
    len = 999;
end

off = 0;

% loop through all trials
for i=1:length(oTrials)
    
    n = oTrials(i).Eye.n;
    
    if oTrials(i).Reward == 1
        
        %align to first frame
        t = oTrials(i).Eye.t(1:n) ...
            - oTrials(i).TrialStartDatapixx ...
            - (oTrials(i).Start(1) - oTrials(i).TrialStart) ;
        
        % filter eye size
        v3 = oTrials(i).Eye.v(3, 1:n);
        v3 = filtfilt(b,a, v3);
        
        %get data at first frame until time after
        n0       = find(t<0, 1, 'last')+1;
        inwindow = n0 : n0+len;
        
        
        % add nan if ther stimuli presented period exceeds the eyelink time
        if n0+len > length(t)
           t = [t nan(1, n0+len - length(t))];
           v3 = [v3 nan(1, n0+len - length(v3))]; 
        end
                
        % get stimuli presented period and add time offset to align to sequence
        t = t(inwindow) + off;
        v3 = v3(inwindow);
        
        if p_flag
            % plot the trajectory
            plot(t, v3, 'Color', col(j,:),...
                'ButtonDownFcn', {@LinePressed, oTrials(i), i});
            hold on;
        end
        
        % assign to trial struct
        oTrials(i).puptraj = v3;
        oTrials(i).window  = j;

        
        %adjust offset
        off = t(end);
        j=j+1;
        
    end
    
    
    if oTrials(i).RewardSize>=0.1 || oTrials(i).Reward == 0 
        j=1;
        off = 0;
    end
    
    if i>1 && ~isempty(strfind(fname, 'all'))
        if oTrials(i).TrialStart - oTrials(i-1).TrialStart > 8  
            j=1;
            off = 0;
        end
        if j==5
            disp('----------------- j == 5 ----------------------')
            fname
            j=1;
        end
    end
    
end

if p_flag
    xlabel('time window of consecutive trials, each aligned to first stimuli presentation');
    ylabel('pupil size'); xlim([0 2]);
    title(fname, 'interpreter','none');
    savefig(h, ['C:\Users\Corinna\Documents\CODE\Analysis\plots\PupilSizeDev\' fname(20:end-4) '_rawppsz.fig']);
    close(h);
end

dat.info = 'only binocular data';


end



function LinePressed(line_h, eventdata, trial, i)
    trial
    i
end






