function [exinfo] = addAll(exinfo)

% % % load list of index
load('Z:\Corinna\SharedCode\Katsu\listpvalue_modulation.mat')
ss = [4 6 8];
for i = 1:length(exinfo)
    
    %% initialization & add p-vals
    exinfo(i).p_modulation = [];        
    try
        exinfo(i).p_modulation = p_modulation{i};
    catch
        disp(['row ' num2str(i) ' is without p_modulation...'])
    end
    
%     %% initialization for the spike isolation analysis
% %     exinfo(i).snr = [];
%     exinfo(i).Lratio = [];
% %     exinfo(i).isolation_score = [];
%     exinfo(i).isolation_distance = [];
% 
%     exinfo(i).isi = [];
%     
% %     exinfo(i).snr_drug = [];
%     exinfo(i).Lratio_drug = [];
% %     exinfo(i).isolation_score_drug = [];
%     exinfo(i).isolation_distance_drug = [];
%     
%     try     
%             [data] = spikeisolation(exinfo(i), 0, 1:14, savefig, savefig);
% %             exinfo(i).snr = data.snr;
%             exinfo(i).isi_fraction = data.fraction;
%             exinfo(i).Lratio = data.Lratio;
% %             exinfo(i).isolation_score = data.isolation_score;
%             exinfo(i).isolation_distance = data.isolation_distance;
%             exinfo(i).isi = [exinfo(i).isi; 1000*data.c1.refp];
%             
%             close all;
%             
%             [data] = spikeisolation(exinfo(i), 1, 1:14, savefig, savefig);
% %             exinfo(i).snr_drug = data.snr;
%             exinfo(i).isi_fraction_drug = data.fraction;
%             exinfo(i).Lratio_drug = data.Lratio;
% %             exinfo(i).isolation_score_drug = data.isolation_score;
%             exinfo(i).isolation_distance_drug = data.isolation_distance;
%             exinfo(i).isi = [exinfo(i).isi; 1000*data.c1.refp];
%             
%             close all;
%             
%             disp(['row ' num2str(i) ' was done with spike isolation analysis (figures are stored) !'])
% 
%     catch
%         disp(['row ' num2str(i) ' was skipped due to error...(spike isolation analysis)'])
%     end
   
    %% variability analysis
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
        exinfo(i).(['ff_re' postfix]) = [];
        exinfo(i).(['nc_re' postfix]) = [];
        exinfo(i).(['sc_re' postfix]) = [];
        exinfo(i).(['sampleSize' postfix]) = [];
        exinfo(i).(['sampleSize_raw' postfix]) = [];
        
        % load data
        try
            ex_su = loadCluster(fname, 'ocul', exinfo(i).ocul, 'loadlfp', false ); % cluster 1 or 2 <=> single unit activity

            fname_mu = strrep(fname, exinfo(i).cluster, 'c0'); % cluster 0 <=> multiunit activity
            ex_mu = loadCluster( fname_mu, 'ocul', exinfo(i).ocul, 'loadlfp', false );
        catch
            disp(['row ' num2str(i) ' is skipped for the variability analysis because data cannot be loaded...'])
            continue
        end

        % sample size used for the variability analysis by CL
        sscl = zeros(1,2);
        stmtype = ex_su.exp.e1.type;
        [u,c] = uniquecount([ex_su.Trials.(stmtype)]);
        sscl(1) = sum(c(1:end-1));
        ex_su2.Trials = getPartialTrials(ex_su.Trials, 2);        
        [~,c] = uniquecount([ex_su2.Trials.(stmtype)]);
        sscl(2) = sum(c(1:end-1));
        exinfo(i).(['sampleSize' postfix]) = [exinfo(i).(['sampleSize' postfix]); sscl];
 
        % subset of trials
        try
            ex_su.Trials(1:20) = [];
            ex_mu.Trials(1:20) = [];
        catch
            disp(['Trials less than 20...'])
            continue
        end
        for s = 1:3
            try                
                [ex_suc, samplesize, c] = getPartialTrialSS(ex_su, ss(s));
                [ex_muc] = getPartialTrialSS(ex_mu, ss(s));                        
            catch
                if d==1 
                    disp(['row ' num2str(i) ' has an error for using the subset with the minimum repeat at ' num2str(ss(s))])
                end
                continue
            end
            if samplesize(2)==0
                 if d==1 
                    disp(['row ' num2str(i) ' has an error for using the subset with the minimum repeat at ' num2str(ss(s))])
                end
                continue
            end

            uc = unique([ex_suc.Trials.(stmtype)]);
            
            % normalize firing data
            [ex_suc, spkrate_su] = znormex(ex_suc, exinfo(i), 0);
            [ex_muc, spkrate_mu] = znormex(ex_muc, exinfo(i), 0);

            % fano factor
            if mean(ismember('all.grating', {exinfo.fname}))==1 && ...
                mean(ismember('all.grating', {exinfo.fname_drug}))==1
                ff = nan(1, length(spkrate_su.cvar));
            else
                ff =  [spkrate_su.cvar]./[spkrate_su.cmn];
            end

            % Noise Correlation
            [rsc, prsc] = corr([ex_suc.Trials.zspkcount]', [ex_muc.Trials.zspkcount]', ...
                'rows', 'complete');

            % Signale Correlation
            [rsig, prsig] = corr([spkrate_su.cmn]', [spkrate_mu.cmn]', ...
                'rows', 'complete');

            % add data
            ffs = nan(1, length(u));
            ffs(ismember(u,uc)) = ff;
            exinfo(i).(['ff_re' postfix]) = [exinfo(i).(['ff_re' postfix]); ffs];
            exinfo(i).(['nc_re' postfix]) = [exinfo(i).(['nc_re' postfix]); [rsc, prsc]];
            exinfo(i).(['sc_re' postfix]) = [exinfo(i).(['sc_re' postfix]); [rsig, prsig]];
            exinfo(i).(['sampleSize' postfix]) = [exinfo(i).(['sampleSize' postfix]); samplesize];
            exinfo(i).(['sampleSize_raw' postfix]) = [exinfo(i).(['sampleSize_raw' postfix]); c];
        end

    end
    
    disp(['row ' num2str(i) ' is done with the variability analysis!'])
end

% exinfo_new = exinfo;
% save('Z:\Katsuhisa\interaction_project\dataset\Data\exinfo.mat', 'exinfo', '-v7.3'); 
% disp('The matrix was saved.')