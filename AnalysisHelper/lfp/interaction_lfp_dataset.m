function interaction_lfp_dataset
% run 'pupil_interaction.m' and 'pupil_lfp.m' in the RC dataset

% load data
if mean(ismember('gpfs0', cd))==1
    load('/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/Data/rcdat.mat') 
else
   load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\rcdat.mat') 
end

pss = struct('session', []);

% loop for pairs of sessions
for i = 1:length(dat)  
    try
        [pss.session(i).psintr_su] = pupil_interaction(dat(i), 'SU');
        pss.session(i).psintr_su_exist = 1;
    catch
        pss.session(i).psintr_su_exist = 0;
        disp(['error in pupil interaction analysis SU at session ' num2str(i)])
    end
    try
        [pss.session(i).psintr_mu] = pupil_interaction(dat(i), 'MU');
        pss.session(i).psintr_mu_exist = 1;
    catch
        pss.session(i).psintr_mu_exist = 0;
        disp(['error in pupil interaction analysis MU at session ' num2str(i)])
    end
    try 
        [pss.session(i).pslfp] = pupil_lfp(dat(i));
        pss.session(i).pslfp_exist = 1;
    catch
        pss.session(i).pslfp_exist = 0;
        disp(['error in pupil vs lfp analysis at session ' num2str(i)])
    end
    disp(['------------------ session ' num2str(i) ' is processed! ----------------------'])
end

% save data structure
if mean(ismember('gpfs0', cd))==1
    save(['/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/Data/pss.mat'], 'pss','-v7.3') 
    disp('output structure pss saved')
else
    save(['Z:\Katsuhisa\serotonin_project\LFP_project\Data\pss.mat'], 'pss','-v7.3') 
    disp('output structure pss saved')
end
