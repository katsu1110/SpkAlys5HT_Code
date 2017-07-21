function col = markerFaceAssignment( exinfo )
% assign face color according to applied drug

if exinfo.pmodulation(2)<0.05
    if exinfo.is5HT
        col = 'r';
    elseif ~exinfo.is5HT
        col = 'k';
    end
else
    col = 'w';
end
    
end

