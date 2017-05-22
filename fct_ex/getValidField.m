function exinfo = getValidField( exinfo )

%---------------------------------------------- choose one trial per neuron
valid  = zeros(1, length(exinfo));
intracmp = zeros(1, length(exinfo));

for neuron = unique([exinfo.id])
    
    neuron_idx = [exinfo.id] == neuron ...
        & ~[exinfo.isadapt];
    
    or_idx = find(strcmp({exinfo.param1}, 'or') & neuron_idx & ~[exinfo.isRC]);
    rc_idx = find(strcmp({exinfo.param1}, 'or') & neuron_idx & [exinfo.isRC]);
    sf_idx = find(strcmp({exinfo.param1}, 'sf') & neuron_idx);
    co_idx = find(strcmp({exinfo.param1}, 'co') & neuron_idx);
    sz_idx = find(strcmp({exinfo.param1}, 'sz') & neuron_idx);
    
    %%%
    if any(or_idx)
        [~, i] = max([exinfo(or_idx).r2reg]);
        valid(or_idx(i)) = 1;
    end
    
    %%%
    if any(rc_idx)
        [~, i] = max([exinfo(rc_idx).r2reg]);
        valid(rc_idx(i)) = 1;
    end
    
    %%%
    if any(sf_idx)
        [~, i] = max([exinfo(sf_idx).r2reg]);
        valid(sf_idx(i)) = 1;
    end
    
    %%%
    if any(co_idx)
        [~, i] = max([exinfo(co_idx).r2reg]);
        valid(co_idx(i)) = 1;
    end
    
    %%%
    if any(sz_idx)
        [~, i] = max([exinfo(sz_idx).r2reg]);
        valid(sz_idx(i)) = 1;
    end
    
    
end


valid        = num2cell(valid==1);
[exinfo.valid] = deal(valid{:});


end

