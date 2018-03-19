function col = getCol4Stim(bool, stim)
% assign color according to stimulus and drug condition

col_5HT = hot(18);
col_NaCl = gray(12);
if bool
    switch stim
        case 'or'
            col = col_5HT(4,:);
        case 'sf'
            col = col_5HT(6,:);
        case 'co'
            col = col_5HT(8,:);
        case 'sz'
            col = col_5HT(10,:);
        case 1
            col = col_5HT(4,:);
        case 2
            col = col_5HT(6,:);
        case 3
            col = col_5HT(8,:);
        case 4
            col = col_5HT(10,:);
    end
else
    switch stim
        case 'or'
            col = col_NaCl(1,:);
        case 'sf'
            col = col_NaCl(3,:);
        case 'co'
            col = col_NaCl(5,:);
        case 'sz'
            col = col_NaCl(7,:);
        case 1
            col = col_NaCl(1,:);
        case 2
            col = col_NaCl(3,:);
        case 3
            col = col_NaCl(5,:);
        case 4
            col = col_NaCl(7,:);
    end
end

end