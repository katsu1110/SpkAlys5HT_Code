function runAlllfp

try
    interaction_lfp_dataset;
catch
    disp('error in interaction_lfp_dataset')
end
try
    fitHMM_dataset;
catch
    disp('error in fitHMM_dataset')
end

