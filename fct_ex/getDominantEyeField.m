function exinfo = getDominantEyeField( exinfo )

%--------------- select dominant eye
isdominant = zeros(1, length(exinfo));

% select the dominant eye, ie the trials that are not binocularly
% stimulated
for neuron = unique([exinfo.idi])
    
    neuron_idx = [exinfo.idi] == neuron & ~[exinfo.isadapt]; %& [exinfo.rsqr4drug] >= expVarThres & ~[exinfo.isRC];
    or_idx = find(strcmp({exinfo.param1}, 'or') & neuron_idx & ~[exinfo.isRC]);
    rc_idx = find(strcmp({exinfo.param1}, 'or') & neuron_idx & [exinfo.isRC]);
    sf_idx = find(strcmp({exinfo.param1}, 'sf') & neuron_idx);
    co_idx = find(strcmp({exinfo.param1}, 'co') & neuron_idx);
    sz_idx = find(strcmp({exinfo.param1}, 'sz') & neuron_idx);
    
    
    if any(or_idx)
        nobino  = find([exinfo(or_idx).ocul] ~= 0);
        [~, uid]  = max(cellfun(@mean, {exinfo(or_idx(nobino)).ratemn}));
        if ~isempty(nobino)
            dom_i = [exinfo.ocul] == exinfo(or_idx(nobino(uid))).ocul ...
                & [exinfo.idi] == neuron & strcmp({exinfo.param1}, 'or');
            isdominant(dom_i) = true;
        end
    end
    
    if any(rc_idx)
        nobino  = find([exinfo(rc_idx).ocul] ~= 0);
        [~, uid]  = max(cellfun(@mean, {exinfo(rc_idx(nobino)).ratemn}));
        if ~isempty(nobino)
            dom_i = [exinfo.ocul] == exinfo(rc_idx(nobino(uid))).ocul ...
                & [exinfo.idi] == neuron & strcmp({exinfo.param1}, 'or');
            isdominant(dom_i) = true;
        end
    end
    
    if any(sf_idx)
        nobino  = find([exinfo(sf_idx).ocul] ~= 0);
        [~, uid]  = max(cellfun(@mean, {exinfo(sf_idx(nobino)).ratemn}));
        if ~isempty(nobino)
            dom_i   = [exinfo.ocul] == exinfo(sf_idx(nobino(uid))).ocul ...
                & [exinfo.idi] == neuron & strcmp({exinfo.param1}, 'sf');
            isdominant(dom_i) = true;
        end
    end
    
    if any(co_idx)
        nobino  = find([exinfo(co_idx).ocul] ~= 0);
        [~, uid]  = max(cellfun(@mean, {exinfo(co_idx(nobino)).ratemn}));
        if ~isempty(nobino)
            dom_i   = [exinfo.ocul] == exinfo(co_idx(nobino(uid))).ocul ...
                & [exinfo.idi] == neuron & strcmp({exinfo.param1}, 'co');
            isdominant(dom_i) = true;
        end
    end
    
    
    if any(sz_idx)
        nobino  = find([exinfo(sz_idx).ocul] ~= 0);
        [~, uid]  = max(cellfun(@mean, {exinfo(sz_idx(nobino)).ratemn}));
        if ~isempty(nobino)
            dom_i   = [exinfo.ocul] == exinfo(sz_idx(nobino(uid))).ocul ...
                & [exinfo.idi] == neuron & strcmp({exinfo.param1}, 'sz');
            isdominant(dom_i) = true;
        end
    end
end


isdominant  = num2cell(isdominant==1);
[exinfo.isdominant] = deal(isdominant{:});


% compare these dominant trials and get the one with the highest firing rate
cmpExp = num2cell(zeros(length(exinfo), 1));
[exinfo.cmpExp] = deal(cmpExp{:});

for uid = unique([exinfo.idi]);
    idx = find([exinfo.idi] == uid);
    [~, maxi] =  max( cellfun(@max, {exinfo(idx).ratemn} ) );
    
    if ~isempty(idx)
        exinfo(idx(maxi)).cmpExp = 1;
    end
end


end

