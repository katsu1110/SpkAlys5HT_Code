function idx = getCritIdx(exinfo, crit)
% returns the index for exinfo that fullfiles the criteria crit in each
% session
% use this file ones for


if strcmp(crit, 'all')
    idx = 1:length(exinfo);
else
    % for 5HT files
    is5HT = find([exinfo.is5HT]);
    idx_5HT = getCritIdx_Helpter(exinfo(is5HT), crit);
    
    
    % for NaCl files
    isNaCl = find(~[exinfo.is5HT]);
    idx_NaCl = getCritIdx_Helpter(exinfo(isNaCl), crit);
    
    % combine indices
    idx = [is5HT(idx_5HT), isNaCl(idx_NaCl)];
end
end


function idx = getCritIdx_Helpter(exinfo, crit)


stim = unique({exinfo.param1});
idx = [];

for stim_i = 1:length(stim)
    for c2 = 0:1
        for rc = 0:1
            switch crit
                case 'best r2'
                    for id = unique([exinfo.id])
                        temp = find([exinfo.id] == id & ...
                            strcmp({exinfo.param1}, stim(stim_i)) & ...
                            [exinfo.isRC] == rc & [exinfo.isc2] == c2);
                        [~, idx_id] = max([exinfo(temp).r2reg]);
                        idx = [idx, temp(idx_id)];
                    end
                case 'first ex'
                    for id = unique([exinfo.id])
                        temp = find([exinfo.id] == id & ...
                            strcmp({exinfo.param1}, stim(stim_i)) & ...
                            [exinfo.isRC] == rc & [exinfo.isc2] == c2);
                        [~, idx_id] = min([exinfo(temp).date]);
                        idx = [idx, temp(idx_id)];
                    end
                case 'highest dose'
                    for id = unique([exinfo.id])
                        temp = find([exinfo.id] == id & ...
                            strcmp({exinfo.param1}, stim(stim_i)) & ...
                            [exinfo.isRC] == rc & [exinfo.isc2] == c2);
                        [~, idx_id] = max([exinfo(temp).dose]);
                        idx = [idx, temp(idx_id)];
                    end
                case 'lowest dose'
                    for id = unique([exinfo.id])
                        temp = find([exinfo.id] == id & ...
                            strcmp({exinfo.param1}, stim(stim_i)) & ...
                            [exinfo.isRC] == rc & [exinfo.isc2] == c2);
                        [~, idx_id] = min([exinfo(temp).dose]);
                        idx = [idx, temp(idx_id)];
                    end
            end
        end
    end
end

end





