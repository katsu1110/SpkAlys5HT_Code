function [exinfo_new] = addvariance(exinfo)

for i = 1:length(exinfo)
    load(exinfo(i).fname)
    data = getEyeVar(ex);
    exinfo(i).variance_eye = data.variance;
    
    load(exinfo(i).fname_drug)
    data = getEyeVar(ex);
    exinfo(i).variance_eye_drug = data.variance;
    
    exinfo(i).variance_eye_1 = [];
    exinfo(i).variance_eye_1_drug = [];
    exinfo(i).variance_eye_2 = [];
    exinfo(i).variance_eye_2_drug = [];
    
    disp(['row ' num2str(i) ' was done!'])
end
exinfo_new = exinfo;

function [fixspan] = getEyeVar(ex)

% initialize output structure ============================================
fixspan = struct('x',[],'y',[],'accuracy',[],'theta',[],'SI',[],'variance',[],'edge',[],'recentering',[],'trial',[]);

% degree per pixels of a display ===============================
% determine which display was used
if isfield(ex,'setup')
        if ex.setup.monitorWidth==56
                screenNum = 1;
        else
                screenNum = 2;                
        end
else
    screenNum = ex.screen_number;
end

% pre-computed Degree Per Pixels
if screenNum==1
%         load('\\172.25.250.112\nienborg_group\Katsuhisa\data\mango\analysis\dpp1.mat')
        dpp = 0.0167;
elseif screenNum==2
%         load('\\172.25.250.112\nienborg_group\Katsuhisa\data\kiwi\analysis\dpp2.mat')
        dpp = 0.0117;
end


% gain of eye calibration =========================================
if isfield(ex.eyeCal,'RXGain') && isfield(ex.eyeCal,'LXGain')
   rgain = [ex.eyeCal.RXGain ex.eyeCal.RYGain];
   lgain = [ex.eyeCal.LXGain ex.eyeCal.LYGain];
   gain = (rgain + lgain)/2;
elseif isfield(ex.eyeCal,'RXGain') && ~isfield(ex.eyeCal,'LXGain')
         gain = [ex.eyeCal.RXGain ex.eyeCal.RYGain];
elseif ~isfield(ex.eyeCal,'RXGain') && isfield(ex.eyeCal,'LXGain')
        gain = [ex.eyeCal.LXGain ex.eyeCal.LYGain];
else
    gain = [ex.eyeCal.XGain ex.eyeCal.YGain];
end


% offset position ===================================================
if isfield(ex.eyeCal,'Delta')
   pos0 = [mean([ex.eyeCal.Delta.RX0]) mean([ex.eyeCal.Delta.RY0])];
elseif isfield(ex.eyeCal,'RX0')
    pos0 = [ex.eyeCal.RX0 ex.eyeCal.RY0];
else
        pos0 = [ex.eyeCal.X0 ex.eyeCal.Y0];
end

% version of ex-file =================================================
if isfield(ex.Trials(1),'TrialStart_remappedGetSecs')
      time = 'N';   % new
elseif ~isfield(ex.Trials(1),'TrialStart_remappedGetSecs')...
        && isfield(ex.Trials(1),'TrialStartDatapixx')...
        && isfield(ex.Trials(1).times, 'lastEyeReadingDatapixx')
      time = 'O';   % old
else
        time = 'C';    % classic
end

% pool eye-position over all the available trials ==============================
avtr = find([ex.Trials.Reward]==1);
len_tr = length(avtr);
Trials = ex.Trials(avtr);
counts = 0;

% convert trial numbers where recentering occured
fixspan.recentering.tr = [];
if isfield(ex.eyeCal,'Delta')
    if isfield(ex.eyeCal.Delta,'TrialNo')
        for i = 1:length([ex.eyeCal.Delta.TrialNo])-1
            if ~isempty(ex.eyeCal.Delta(i).TrialNo)
                delta = ex.eyeCal.Delta(i).TrialNo - avtr;
                fixspan.recentering.tr = [fixspan.recentering.tr ex.eyeCal.Delta(i).TrialNo - min(delta(delta > 0))];
            end
        end
    end
end

% initialization
start = 1;
eyepos.trial = struct('x', [], 'y', []);
fixspan.recentering.tr = unique([0 fixspan.recentering.tr]);
fixspan.recentering.window_offsetX = dpp*pos0(1);
fixspan.recentering.window_offsetY = dpp*pos0(2);
fixspan.recentering.eye_offsetX = [];
fixspan.recentering.eye_offsetY = [];

for i = 1:len_tr
        
        % timing information in a trial (start and stop of the stimulus)
        if time == 'N'
            t = Trials(i).Eye.t(1:Trials(i).Eye.n)-Trials(i).TrialStartDatapixx;
            st = Trials(i).Start - Trials(i).TrialStart_remappedGetSecs;            

            [~,stpos] = min(abs(t-st(1)));
            [~,enpos] = min(abs(t-st(end))); 
              
        elseif time == 'O'
            delta = Trials(i).Eye.t(Trials(i).Eye.n) - Trials(i).TrialStartDatapixx - Trials(i).times.lastEyeReading;
            t = Trials(i).Eye.t(1:Trials(i).Eye.n)-Trials(i).TrialStartDatapixx-delta;
            st = Trials(i).Start - Trials(i).TrialStart;

            [~,stpos] = min(abs(t-st(1)));
            [~,enpos] = min(abs(t-st(end)));
            
         elseif time == 'C'
             stpos = floor(Trials(i).times.startFixation*sampRate);
             enpos = floor((Trials(i).Start(end) - Trials(i).Start(1))*sampRate);
        end

        % raw x & y eye-positions
        if length(Trials(i).Eye.v(:,1)) > 3     % binocular eye-tracking
                xeye = mean(Trials(i).Eye.v([1 4], stpos:enpos),1);        
                yeye = mean(Trials(i).Eye.v([2 5], stpos:enpos),1);        
        else     % monocular eye-tracking
                xeye = mean(Trials(i).Eye.v(1, stpos:enpos),1);  
                yeye = mean(Trials(i).Eye.v(2, stpos:enpos),1);  
        end
        
        % convert voltage into vda      
        xeye = (xeye - pos0(1))*gain(1)*dpp;
        yeye = (yeye - pos0(2))*gain(2)*dpp; 
%         xeye = (xeye - xpos0)*gain(1)*dpp;
%         yeye = (yeye - ypos0)*gain(2)*dpp;   

        % adjust centering?
        if isfield(ex.eyeCal,'Delta')
            if isfield(ex.eyeCal.Delta,'TrialNo')
               if ismember(i ,fixspan.recentering.tr)
                   counts = counts + 1;                  
                   fixspan.recentering.window_offsetX = ...
                           [fixspan.recentering.window_offsetX dpp*ex.eyeCal.Delta(counts).RX0];
                   fixspan.recentering.window_offsetY = ...
                           [fixspan.recentering.window_offsetY dpp*ex.eyeCal.Delta(counts).RY0];
                   fixspan.recentering.eye_offsetX = ...
                           [fixspan.recentering.eye_offsetX median(fixspan.x(start:end))];
                   fixspan.recentering.eye_offsetY = ...
                           [fixspan.recentering.eye_offsetY median(fixspan.y(start:end))];
                   start = length(fixspan.x) + 1;
               end
            end
        end        
           
        % pool the eye-position across trials    
        eyepos.trial(i).x = xeye;
        eyepos.trial(i).y = yeye;
        if length(Trials(i).Eye.v(:,1)) > 3     % binocular eye-tracking
                fixspan.x = [fixspan.x, xeye];        
                fixspan.y = [fixspan.y, yeye];        
        else     % monocular eye-tracking
                fixspan.x = [fixspan.x, xeye];  
                fixspan.y = [fixspan.y, yeye];  
        end
    
end

% equalizing the length of 'window_offset' and 'eye_offset'
fixspan.recentering.eye_offsetX = ...
        [fixspan.recentering.eye_offsetX median(fixspan.x(start:end))];
fixspan.recentering.eye_offsetY = ...
        [fixspan.recentering.eye_offsetY median(fixspan.y(start:end))];

% remove offset of eye-positions
fixspan.x = fixspan.x - median(fixspan.x);
fixspan.y = fixspan.y - median(fixspan.y);

% structualize =======================================
fixspan.variance = [var(fixspan.x) var(fixspan.y)];

