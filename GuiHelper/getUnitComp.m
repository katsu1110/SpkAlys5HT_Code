function exinfo_out = getUnitComp(spec, exinfo)
% reduces the exinfo structure to those entries who are relevant and are in
% accordance to the GUI call. These are the ones complying to the input
% structure spec. The structure contains information about the desired
% ocular condition and desired stimulus experiments (or, co, ...).
% it also ensured that each unit is represented only once per drug
% condition.
% 
% 
% @CL



% In an older version, it was possible to plot the dominant versus the
% non-dominant eye condition or between one and another stimulus condition. 
% This is checked in the first to if-cases.
% I did not use this for a long time so it could be that it is not working
% anymore.
if strcmp(spec.stimx, 'all stimuli cond') && ~(strcmp(spec.stimy, spec.stimx))
    ix = findCorrespStimIdx(spec.stimy, spec.eyex, exinfo, 0);
    iy = findCorrespStimIdx(spec.stimy, spec.eyey, exinfo, 1);
elseif strcmp(spec.stimy, 'all stimuli cond') && ~(strcmp(spec.stimy, spec.stimx))
    ix = findCorrespStimIdx(spec.stimx, spec.eyex, exinfo, 1);
    iy = findCorrespStimIdx(spec.stimx, spec.eyey, exinfo, 0);
else
    % comparing across experiments with the same stimulus variation
    ix = findCorrespStimIdx(spec.stimx, spec.eyex, exinfo);
    iy = findCorrespStimIdx(spec.stimy, spec.eyey, exinfo);
end


if all(ix == iy)
    exinfo_out = singleUnitsOnly(exinfo(ix | iy));
else 
    error('the indices are not matching. check the stimulus condition');
end


end



function ind = findCorrespStimIdx(specstim, speceye, exinfo, neg_flag)
% find the stimulus and ocular specific entries in exinfo
%
% specstim is the experiment specification (or, co, RC, adapt, ...)
% speceye is the ocular condition 
% exinfo is the result structure
% neg_flag indicates whether the c
%
%
%



if nargin < 4
    neg_flag = false;
end

% ocular condition
switch speceye
    case 'all'
        eye_spec = zeros(1, length(exinfo));
        for idi = unique([exinfo.idi])
            
            % find the ocular condition that elicited the highest response
            idx = find([exinfo.idi] == idi);
            [~, maxi] =  max( cellfun(@max, {exinfo(idx).ratemn} ) );
            
            if ~isempty(idx)
                eye_spec(idx(maxi)) = 1;
            end
        end
        
    case 'dominant eye'
        eye_spec = [exinfo.isdominant];
    case 'non-dominant eye'
        eye_spec = ~[exinfo.isdominant] & [exinfo.ocul]~=0;
end


if neg_flag
    if strcmp(specstim , 'RC')
        ind = ~[exinfo.isRC] & eye_spec;
    elseif strcmp(specstim , 'adapt')
        ind = ~[exinfo.isadapt] & eye_spec;
    else
        ind = ~strcmp({exinfo.param1}, specstim) & eye_spec ...
            & ~[exinfo.isadapt] & ~[exinfo.isRC]; 
    end
else
    if strcmp(specstim , 'all stimuli cond')
        ind = eye_spec;
    elseif strcmp(specstim , 'RC')
        ind = [exinfo.isRC];
    elseif strcmp(specstim , 'adapt')
        ind = [exinfo.isadapt] & eye_spec;
    else
        ind = strcmp({exinfo.param1}, specstim) & eye_spec ...
            & ~[exinfo.isadapt] & ~[exinfo.isRC]; 
    end
end



end



function expInfo_out = singleUnitsOnly(exinfo)
% if there are conflicting data, i.e. experiments from the same unit, the
% one with highest baseline response is used


% make sure to distinguish between clusters and different drug applications
idx_5HT1 = singleUnitsOnlyHelper(exinfo, [exinfo.is5HT]& [exinfo.isc2]);
idx_5HT2 = singleUnitsOnlyHelper(exinfo, [exinfo.is5HT]& ~[exinfo.isc2]);
idx_NaCl1 = singleUnitsOnlyHelper(exinfo, ~[exinfo.is5HT]& [exinfo.isc2]);
idx_NaCl2 = singleUnitsOnlyHelper(exinfo, ~[exinfo.is5HT]& ~[exinfo.isc2]);

expInfo_out = exinfo([idx_5HT1; idx_5HT2; idx_NaCl1; idx_NaCl2]);
    

end


function idx_su = singleUnitsOnlyHelper(exinfo, id_5HT)
% find the unit with the highest baseline response

idx_su = [];
id_i = unique([exinfo(id_5HT).id]);
for i = 1:length(id_i)
    
    idx_temp = find([exinfo.id] == id_i(i) & id_5HT); % indices of current unit
    
    if length(idx_temp) > 1
        [~,k] = max(cellfun(@max, {exinfo(idx_temp).ratemn}));
        idx_su = [idx_su; idx_temp(k)];
    else
        idx_su = [idx_su; idx_temp];
    end
end

end