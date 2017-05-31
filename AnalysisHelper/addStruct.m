function exinfo = addStruct( exinfo )
% post hoc expansion of exinfo 
%
% copy the fitting quality to the field exinfo.gaussr2.
% 
% @CL



for i = 1:length(exinfo)
    
    if ~isempty(exinfo(i).fitparam)
        if isfield(exinfo(i).fitparam, 'OR') && isfield(exinfo(i).fitparam, 'r2')
            
            maxi = find( max(exinfo(i).sdfs.y(1,:)) == max(exinfo(i).sdfs.y(1,:)), 1, 'first');
            exinfo(i).gaussr2 = exinfo(i).fitparam.OR(maxi).r2; 
            
        elseif isfield(exinfo(i).fitparam, 'r2')
            exinfo(i).gaussr2 = exinfo(i).fitparam.r2;
        else
            exinfo(i).gaussr2 = 0;
        end
    else
        exinfo(i).gaussr2 = 0;
    end
       
    
    
    if ~isempty(exinfo(i).fitparam_drug)
        if isfield(exinfo(i).fitparam_drug, 'OR') && isfield(exinfo(i).fitparam_drug, 'r2')
            
            maxi = find( max(exinfo(i).sdfs_drug.y(1,:)) == max(exinfo(i).sdfs_drug.y(1,:)), 1, 'first');
            exinfo(i).gaussr2 = exinfo(i).fitparam_drug.OR(maxi).r2;
            
        elseif isfield(exinfo(i).fitparam_drug, 'r2')
            exinfo(i).gaussr2_drug = exinfo(i).fitparam_drug.r2;
        else
            exinfo(i).gaussr2_drug = 0;
        end
    else
        exinfo(i).gaussr2_drug = 0;
    end
    
    
    if isnan(exinfo(i).gaussr2)
        exinfo(i).gaussr2 = 0;
    end
    
    if isnan(exinfo(i).gaussr2_drug)
        exinfo(i).gaussr2_drug = 0;
    end
    
  

end


end

