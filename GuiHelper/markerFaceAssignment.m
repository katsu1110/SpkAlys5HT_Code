function col = markerFaceAssignment( exinfo )
% assign face color according to applied drug

if isfield(exinfo, 'p_modulation')
%     disp(['p_val: ' num2str(exinfo.p_modulation{1}(2))])
    if exinfo.p_modulation(2) >= 0.05
        col = 'w';
    else
        if exinfo.is5HT
            col = 'r';
        elseif ~exinfo.is5HT
            col = 'k';
        end
    end
else
    if exinfo.is5HT
            col = 'r';
        elseif ~exinfo.is5HT
            col = 'k';
    end
end

