function col = markerFaceAssignment( exinfo )
% assign face color according to applied drug

if exinfo.is5HT
    col = 'r';
elseif ~exinfo.is5HT
    col = 'k';
end
    
end

