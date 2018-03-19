function [LFPinfo] = addLFP(exinfo, pupil_flag)
%% run 'LFPanalyzer.m' in batch using exinfo
%
%
% written by Katsuhisa (29.09.17)
% ++++++++++++++++++++++++++++++++

LFPinfo = struct('session', []);
for i = 1:length(exinfo)
    try
        if pupil_flag==1
            [para] = LFPanalyzer(exinfo(i),'plot','pupil','save');
        else
            [para] = LFPanalyzer(exinfo(i),'plot','save');
        end
        LFPinfo.session(i).exist = 1;
    catch
        LFPinfo.session(i).exist = 0;
    end
    if LFPinfo.session(i).exist == 1
        LFPinfo.session(i).results = para;
        disp(['row ' num2str(i) ' was analyzed and figures were stored.'])
    else
        disp(['row ' num2str(i) ' had an error.'])
    end
end
    
if pupil_flag==1
    save('Z:\Katsuhisa\LFP_project\Data\LFPinfo_pupil.mat', 'LFPinfo', '-v7.3')
else
    save('Z:\Katsuhisa\LFP_project\Data\LFPinfo.mat', 'LFPinfo', '-v7.3')
end
disp('LFPinfo was saved!')