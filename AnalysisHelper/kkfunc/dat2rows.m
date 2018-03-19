function [rows] = dat2rows(dat, exinfo)

len_dat = length(dat.x);
len_exinfo = length(exinfo);
rows = nan(1, len_dat);

for i = 1:len_dat
    for e = 1:len_exinfo
        if strfind(dat.expInfo(i).fname_drug, exinfo(e).fname_drug)
            rows(i) = e;
            break
        end
    end
end

