function exinfo = getValidField( exinfo )
% exinfo = getValidField( exinfo )
% 
% adds a new field valid that identifies the experiment pair with the best
% regression fit for each unit across the same stimulus dimension.
% 
% 
% @CL
%


%intitialize the variable carrying the boolean values that identify the
%experiment pairs with the best regression fits per unit and tested
%stimulus dimension
hasbestregfit  = false(1, length(exinfo));

% loop through all units
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
        hasbestregfit(or_idx(i)) = 1;
    end
    
    %%%
    if any(rc_idx)
        [~, i] = max([exinfo(rc_idx).r2reg]);
        hasbestregfit(rc_idx(i)) = 1;
    end
    
    %%%
    if any(sf_idx)
        [~, i] = max([exinfo(sf_idx).r2reg]);
        hasbestregfit(sf_idx(i)) = 1;
    end
    
    %%%
    if any(co_idx)
        [~, i] = max([exinfo(co_idx).r2reg]);
        hasbestregfit(co_idx(i)) = 1;
    end
    
    %%%
    if any(sz_idx)
        [~, i] = max([exinfo(sz_idx).r2reg]);
        hasbestregfit(sz_idx(i)) = 1;
    end
    
    
end

% assign the best regression fit identifiers as field to the result
% structure
hasbestregfit        = num2cell(hasbestregfit==1);
[exinfo.valid] = deal(hasbestregfit{:});


end

