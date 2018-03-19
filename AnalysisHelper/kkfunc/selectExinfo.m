function [exinfo0, exinfo1] = selectExinfo(exinfo)


rank_thre = 2;
iso_thre = 20;
L_thre = 0.05;
len_exp = length(exinfo);
out0 = zeros(1, len_exp);
out1 = zeros(1, len_exp);
for i = 1:len_exp
    rank = min([exinfo(i).spkqual_base, exinfo(i).spkqual_drug]);
    if rank > rank_thre
        out0(i) = 1;
    end
    
    try 
        isodist = min([exinfo(i).isolation_distance(1), exinfo(i).isolation_distance_drug(1)]);
        lratio = max([exinfo(i).Lratio(1), exinfo(i).Lratio_drug(1)]);
        if isodist < iso_thre || lratio >= L_thre
            out1(i) = 1;
        end
    catch
        out1(i) = 1;
    end
end

exinfo0 = exinfo(out0==0);
exinfo1 = exinfo(out1==0);