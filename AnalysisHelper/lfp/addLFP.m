function [LFPinfo] = addLFP(exinfo)
%% run 'LFPanalyzer.m' in batch using exinfo
%
% load('Z:\Katsuhisa\serotonin_project\dataset\Data\exinfo.mat')
% written by Katsuhisa (29.09.17)
% ++++++++++++++++++++++++++++++++

load('Z:\Corinna\SharedCode\Katsu\list_RC.mat')
load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stim_cond_2007.mat')
LFPinfo = struct('session', []);
for i = 1:length(exinfo)
    % good unit or not
    if (exinfo(i).isRC==1 && ismember(i, list_RC)) || ...
         (exinfo(i).isRC==0 && ismember(i, incl_i))
        LFPinfo.session(i).goodunit = 1;
    else
        LFPinfo.session(i).goodunit = 0;
    end
    % LFP analysis if the session's stimulus is either 'co' or 'RC'
    if strcmp(exinfo(i).param1, 'co') || exinfo(i).isRC==1
        try
            [para] = LFPanalyzer(exinfo(i),'plot','save');
            LFPinfo.session(i).exist = 1;
        catch
            LFPinfo.session(i).exist = 0;
        end
        if LFPinfo.session(i).exist == 1
            LFPinfo.session(i).results = para;
            disp(['Session ' num2str(i) ' was analyzed and figures were stored.'])
        else
            disp(['Session ' num2str(i) ' had an error.'])
        end
    else
        LFPinfo.session(i).exist = 0;
        disp(['Session ' num2str(i) ' is neither CO nor RC experiment. Skipped...'])
    end
end
    
save('Z:\Katsuhisa\serotonin_project\LFP_project\Data\LFPinfo.mat', 'LFPinfo', '-v7.3')
disp('LFPinfo was saved!')