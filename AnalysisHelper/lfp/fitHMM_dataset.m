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

binsize = 10;
% binsize = 15;
hmms = struct('session', []);
n_comp = [1,2,3];

% loop for pairs of sessions
for b = 1:length(binsize)
    disp(['bin size: ' num2str(binsize(b))])
    for s = 1:2 % SU or MU
        if s==1
            disp(['-------- single unit ----------'])
        else
            disp(['-------- multi unit ----------'])
        end
        for i = 1:length(dat)
            % files
            if s==1 % single-unit
                fname = dat(i).fname;
                fname_drug = dat(i).fname_drug;
            else % multi-unit
                fname = strrep(dat(i).fname, dat(i).cluster, 'c0'); % cluster 0 <=> multiunit activity
                fname_drug = strrep(dat(i).fname_drug, dat(i).cluster, 'c0'); % cluster 0 <=> multiunit activity
            end
            
            % control            
            ex = loadCluster(fname, 'ocul', dat(i).ocul, 'loadlfp', false);
            spikecount = ex2sc(ex, binsize(b));        

            % drug           
            ex = loadCluster(fname_drug, 'ocul', dat(i).ocul, 'loadlfp', false);
            spikecount_drug = ex2sc(ex, binsize(b));

            % into structure
            hmms.session(i).unit(s).fname = fname;
            hmms.session(i).unit(s).fname_drug = fname_drug;
            hmms.session(i).unit(s).binsize = binsize(b);
            hmms.session(i).unit(s).spikecount = spikecount;
            hmms.session(i).unit(s).spikecount_drug = spikecount_drug;

            % fit GPFA
            try
                hmms.session(i).unit(s).gpfa = fitGPFA(spikecount);
                disp(['GPFA fitted on session : ' num2str(i)])
            catch
                disp(['error in GPFA fit on session : ' num2str(i)])
            end
            try
                hmms.session(i).unit(s).gpfa_drug = fitGPFA(spikecount_drug);
                disp(['GPFA fitted on session : ' num2str(i)])
            catch
                disp(['error in GPFA fit on session : ' num2str(i)])
            end

            % fit HMM
            for j = 1:n_comp(end)
                try
                    [hmms.session(i).unit(s).estimate(j).fit] = fitHMM(spikecount, n_comp(j));
                    hmms.session(i).unit(s).estimate(j).exist = 1;
                    disp(['HMM fitted on session : ' num2str(i)  ' control data: N =  ' num2str(j)])
                catch
                    hmms.session(i).unit(s).estimate(j).exist = 0;
                    disp(['HMM fit; error on session : ' num2str(i)  ' control data: N =  ' num2str(j)])
                end
                try
                    [hmms.session(i).unit(s).estimate_drug(j).fit] = fitHMM(spikecount_drug, n_comp(j));
                    hmms.session(i).unit(s).estimate_drug(j).exist = 1;
                    disp(['HMM fitted on session : ' num2str(i)  ' drug data: N =  ' num2str(j)])
                catch                
                    hmms.session(i).unit(s).estimate_drug(j).exist = 0;
                    disp(['HMM fit; error on session : ' num2str(i)  ' drug data: N =  ' num2str(j)])
                end
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
% 350ms after the stimulus onset
start = 0.35;
ncol = round((1000/b)*(ex.fix.duration - start));
sc = zeros(len_tr, ncol);
for i = 1:len_tr    
    begin = start;
    for c = 1:ncol
        sc(i,c) = sum(ex.Trials(i).Spikes > begin & ex.Trials(i).Spikes <= begin + 0.001*b);
        begin = begin + 0.001*b;
    end
end
