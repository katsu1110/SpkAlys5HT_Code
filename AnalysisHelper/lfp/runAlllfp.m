function runAlllfp
 
% check_FRonLFP
try
    interaction_lfp_dataset;
catch
    disp('error in interaction_lfp_dataset')
end
try
    if mean(ismember('gpfs0', cd))==1
        load('/gpfs01/nienborg/group/Katsuhisa/serotonin_project/dataset/Data/exinfo.mat')
    else
        load('Z:/Katsuhisa/serotonin_project/dataset/Data/exinfo.mat')
    end
    addLFP(exinfo);
catch
    disp('error in addLFP')
end
% try
%     fitHMM_dataset;
% catch
%     disp('error in fitHMM_dataset')
% end
