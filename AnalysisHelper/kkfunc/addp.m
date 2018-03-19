function [exinfo_cl] = addp(exinfo)

 load('Z:\Corinna\SharedCode\Katsu\listpvalue_modulation.mat')
for i = 1:length(exinfo)
        
        % initialization
        exinfo(i).p_modulation = [];
        
        % load data
        try
            exinfo(i).p_modulation = p_modulation{i};

        catch
            disp(['row ' num2str(i) ' is skipped because of errors...'])
            continue
        end
   
    disp(['row ' num2str(i) ' is done!'])
end
exinfo_cl = exinfo;
% save('Z:\Katsuhisa\interaction_project\dataset\Data\exinfo_cl.mat', 'exinfo_cl', '-v7.3'); 

