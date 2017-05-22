
function oTrials = addWindow(exinfo, oTrials)

wind = 1;

% loop through all trials and add an additional field window that gives the
% number in the 2s stimuli sequence 
for i=1:length(oTrials)
        
    if oTrials(i).Reward == 1
        oTrials(i).window = wind;
        wind=wind+1;
    else
        oTrials(i).window = 0;
        wind=1;
    end
    
    
    % if the trial sequence was rewarded OR this trial showed a blank,
    % reset window index wind
    if oTrials(i).RewardSize>=0.1 || oTrials(i).(exinfo.param1)>1000
        wind=1;
    end
    
    if i>1 && ~isempty(strfind(exinfo.fname, 'all'))
        if oTrials(i).TrialStart - oTrials(i-1).TrialStart > 8  || wind>=5
            wind=1;
        end
    end
end

end