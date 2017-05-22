
function fname = voltList(exinfo)

for i = 1:length(exinfo)
    fprintf('session %1.1f \n', exinfo(i).id);
    fname{i} = strrep(exinfo(i).fname_drug, 'D', 'Z');
    try
        load(fname{i});
        volt(i) = getVolt(ex);
    catch
        volt(i) = 1;
    end
    
end

fname = fname(isnan(volt));
fname = unique(fname');
end


