function [isodist_base, isodist_drug] = isolationDistVec(exinfo, type)

switch type 
    case 1
        name = 'isolation_distance';
    case 2
        name = 'Lratio';
end

len_ex = length(exinfo);
isodist_base = nan(1, len_ex);
isodist_drug = nan(1, len_ex);
for i = 1:len_ex
    if ~isempty(exinfo(i).(name))
        isodist_base(i) = exinfo(i).(name)(2);
    end
    if ~isempty(exinfo(i).([name '_drug']))
        isodist_drug(i) = exinfo(i).([name '_drug'])(2);
    end
end

isodist_base(isnan(isodist_base)) = nanmedian(isodist_base);
isodist_drug(isnan(isodist_drug)) = nanmedian(isodist_drug);
