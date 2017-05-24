function volt = getVolt(ex)
% volt = getVolt(ex)
% returns the recorded voltage from the Header.
% When ex contains concatenated experiments, the first experiment is
% considered.
% 
% @CL


    
volt = nan; % default  

if isfield(ex.Header, 'iontophoresisVoltage')
    volt = ex.Header.iontophoresisVoltage;
    
elseif isfield(ex.Header, 'Headers')
    if isfield(ex.Header.Headers(1), 'iontophoresisVoltage')
        volt = ex.Header.Headers(1).iontophoresisVoltage;
    end
    
end

end