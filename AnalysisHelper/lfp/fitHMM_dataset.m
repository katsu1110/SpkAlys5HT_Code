function [hmms] = fitHMM_dataset(session_idx)
% fit HMM for the RC datasets
% INPUT: session_idx ... vector of integers: 1 - 44 or 'all'

% load data
if mean(ismember('gpfs0', cd))==1
    load('/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/Data/rcdat.mat') 
else
   load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\rcdat.mat') 
end

if nargin < 1 || strcmp(session_idx, 'all')
    session_idx = 1:length(dat);
end
dat = dat(session_idx);

binsize = 10:10:50;
% binsize = 15;
hmms = struct('session', []);
n_comp = [1,2];

% loop for pairs of sessions
for b = 1:length(binsize)
    disp(['bin size: ' num2str(binsize(b))])
    for i = 1:length(dat)  
        % files
        fname_mu = strrep(dat(i).fname, dat(i).cluster, 'c0'); % cluster 0 <=> multiunit activity
        fname_mu_drug = strrep(dat(i).fname_drug, dat(i).cluster, 'c0'); % cluster 0 <=> multiunit activity
        
        % control            
        ex_mu = loadCluster(fname_mu, 'ocul', dat(i).ocul, 'loadlfp', false);
        spikecount = ex2sc(ex_mu, binsize(b));        

        % drug           
        ex_mu = loadCluster(fname_mu_drug, 'ocul', dat(i).ocul, 'loadlfp', false);
        spikecount_drug = ex2sc(ex_mu, binsize(b));

        % into structure
        hmms.session(i).fname = fname_mu;
        hmms.session(i).fname_drug = fname_mu_drug;
        hmms.session(i).binsize = binsize(b);
        hmms.session(i).spikecount = spikecount;
        hmms.session(i).spikecount_drug = spikecount_drug;
        
        % fit HMM
        for j = 1:n_comp(end)
            try
                [hmms.session(i).estimate(j).fit] = fitHMM(spikecount, n_comp(j));
                hmms.session(i).estimate(j).exist = 1;
                disp(['HMM fitted on session : ' num2str(i)  ' control data: N =  ' num2str(j)])
            catch
                hmms.session(i).estimate(j).exist = 0;
                disp(['error on session : ' num2str(i)  ' control data: N =  ' num2str(j)])
            end
            try
                [hmms.session(i).estimate_drug(j).fit] = fitHMM(spikecount_drug, n_comp(j));
                hmms.session(i).estimate_drug(j).exist = 1;
                disp(['HMM fitted on session : ' num2str(i)  ' drug data: N =  ' num2str(j)])
            catch                
                hmms.session(i).estimate_drug(j).exist = 0;
                disp(['error on session : ' num2str(i)  ' drug data: N =  ' num2str(j)])
            end
        end        
    end

    % save data structure
    if mean(ismember('gpfs0', cd))==1
        save(['/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/Data/hmms' num2str(binsize(b)) '.mat'], 'hmms','-v7.3') 
    else
        save(['Z:\Katsuhisa\serotonin_project\LFP_project\Data\hmms' num2str(binsize(b)) '.mat'], 'hmms','-v7.3') 
    end
end

function sc = ex2sc(ex, b)
% spike counts
len_tr = length(ex.Trials);
ncol = round((1000/b)*ex.fix.duration);
sc = zeros(len_tr, ncol);
for i = 1:len_tr
    begin = 0;
    for c = 1:ncol
        sc(i,c) = sum(ex.Trials(i).oSpikes > begin & ex.Trials(i).oSpikes <= begin + 0.001*b);
        begin = begin + 0.001*b;
    end
end