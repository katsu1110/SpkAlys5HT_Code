function [exinfo_var] = addSS(exinfo)


ss = [4 6 8 10];
for i = 1:length(exinfo)
    for d = 1:2
        switch d
            case 1
                fname = exinfo(i).fname;
                postfix = '';                
            case 2
                fname = exinfo(i).fname_drug;
                postfix = '_drug';
        end
        
        % initialization
        exinfo(i).(['sampleSize' postfix]) = [];
        
        % load data
        try
            ex_su = loadCluster(fname, 'ocul', exinfo(i).ocul, 'loadlfp', false ); % cluster 1 or 2 <=> single unit activity

        catch
            disp(['row ' num2str(i) ' is skipped because data cannot be loaded...'])
            continue
        end

        % sample size used for the variability analysis by CL
        sscl = zeros(1,2);
        stmtype = ex_su.exp.e1.type;
        [~,c] = uniquecount([ex_su.Trials.(stmtype)]);
        sscl(1) = sum(c(1:end-1));
        ex_su2.Trials = getPartialTrials(ex_su.Trials, 2);        
        [~,c] = uniquecount([ex_su2.Trials.(stmtype)]);
        sscl(2) = sum(c(1:end-1));
        exinfo(i).(['sampleSize' postfix]) = [exinfo(i).(['sampleSize' postfix]); sscl];

        % subset of trials
        try
            ex_su.Trials(1:20) = [];
        catch
            disp(['Trials less than 20...'])
            continue
        end
        for s = 1:4
            try                
                [~, samplesize] = getPartialTrialSS(ex_su, ss(s));
            catch
                if d==1
                    disp(['row ' num2str(i) ' has an error for using the subset with the minimum repeat at ' num2str(ss(s))])
                end
                continue
            end
            
            % add data
            exinfo(i).(['sampleSize' postfix]) = [exinfo(i).(['sampleSize' postfix]); samplesize];
        end

    end
    
    disp(['row ' num2str(i) ' is done!'])
end
exinfo_var = exinfo;
save('Z:\Katsuhisa\interaction_project\dataset\Data\exinfo_var.mat', 'exinfo_var', '-v7.3'); 

