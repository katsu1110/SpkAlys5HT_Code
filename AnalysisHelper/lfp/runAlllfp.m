function runAlllfp
 
check_FRonLFP
try
    interaction_lfp_dataset;
catch
    disp('error in interaction_lfp_dataset')
end
try
    load('Z:\Katsuhisa\serotonin_project\dataset\Data\exinfo.mat')
    addLFP(exinfo);
catch
    disp('error in addLFP')
end
% try
%     fitHMM_dataset;
% catch
%     disp('error in fitHMM_dataset')
% end

