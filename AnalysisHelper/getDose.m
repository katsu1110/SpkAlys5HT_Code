function [dose, dosernd] = getDose(ex)
%[dose, dosebin] = getDose(ex)
% 
% retract the dose, i.e. the applied current, from the ex Header
% Additionally return the binned version given by getDosebin(dose).
% 
% 
% @CL

dose = -1; % default

if isfield(ex.Header, 'Headers') 
    if isfield(ex.Header.Headers(1), 'iontophoresisEjectionCurrent')
        dose = ex.Header.Headers(1).iontophoresisEjectionCurrent;
    end
else
    if isfield(ex.Header, 'iontophoresisEjectionCurrent')
        dose = ex.Header.iontophoresisEjectionCurrent;
    end
end

if nargout == 2
    dosernd = getDosebin(dose(1));
end

end


%%
function doseend = getDosebin(dose)
% arbitrary binning of applied drug dose

if dose<0
    doseend=-1;
elseif dose==0.1
    doseend=-1;
elseif dose<=5
    doseend=5;
elseif dose<=10
    doseend=10;
elseif dose<=16
    doseend=16;
elseif dose<=25
    doseend=25;
elseif dose<=35
    doseend = 35;
elseif dose<=inf
    doseend = 45;
else
    doseend = dose;
end

end