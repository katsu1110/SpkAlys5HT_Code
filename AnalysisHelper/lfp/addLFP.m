function [LFPinfo] = addLFP(exinfo)
%% run 'LFPanalyzer.m' in batch using exinfo
%
% load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\dataset\Data\exinfo.mat')
% written by Katsuhisa (29.09.17)
% ++++++++++++++++++++++++++++++++

LFPinfo = struct('session', []);
for i = 1:length(exinfo)
    if strcmp(exinfo(i).params(1), 'co') || exinfo(i).isRC==1
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
        disp(['Session ' num2str(i) ' is neither CO nor RC experiment. Skipped...'])
    end
end
    
save('Z:\Katsuhisa\LFP_project\Data\LFPinfo.mat', 'LFPinfo', '-v7.3')
disp('LFPinfo was saved!')