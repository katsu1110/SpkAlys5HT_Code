function c = getCol(ex)

if isfield(ex, 'is5HT')
    if ex.is5HT
        c = 'r';
    else
        c = 'k';
    end
else
    
        
    fname = getFname(ex);
    
    if isempty(strfind(fname, '5HT'))
        c = 'k';
    else
        c = 'r';
    end
    
    
end

end



