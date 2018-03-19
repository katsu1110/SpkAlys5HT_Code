len_ex = length(exinfo);

%% look for serotonin experiment with recovery
list_ser_recover = [];
dat5ht = find([dat.is5HT]==1);

for e = 1:length(exinfo)
    for i = 1:length(dat5ht)
        if isequal(exinfo(e).fname_drug, dat.expInfo(dat5ht(i)).fname_drug)
            list_ser_recover = [list_ser_recover, e];
        end
    end
end

%% look for NaCl
load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stim_cond_2007.mat')
list_nacl = [];
for i = 1:length(incl_i)
    if strcmp(exinfo(incl_i(i)).drugname, 'NaCl')
        list_nacl = [list_nacl, incl_i(i)];
    end
end

%% look for serotonin
list_ser = incl_i(~ismember(incl_i, list_nacl));
list_ser_alone = [];
for i = 1:length(list_ser)
    c = 0;
    for e = 1:len_ex
        if exinfo(e).id==exinfo(list_ser(i)).id
            c = c + 1;
        end
    end
    if c == 1
        list_ser_alone = [list_ser_alone, list_ser(i)];
    end
end


unit_ser_recover = [dat.expInfo([dat.is5HT]==1).id];
unit_ser = [exinfo(list_ser).id];
list_nacl_with = [];
list_nacl_alone = [];
for i = 1:length(list_nacl)
    if ismember(exinfo(list_nacl(i)).id, unit_ser_recover)
        list_nacl_with = [list_nacl_with, list_nacl(i)];
    end
    if ~ismember(exinfo(list_nacl(i)).id, unit_ser)
        list_nacl_alone = [list_nacl_alone, list_nacl(i)];
    end
end


% a) for which no recovery experiment was recorded
% b) for which we had a recovery experiment but the unit didn't achieve full recovery.
nacl_recov_exp = nan(1, length(list_nacl));
nacl_recovered = nan(1, length(list_nacl));
for i = 1:length(list_nacl)
    if exinfo(list_nacl(i)).id==174.5
        nacl_recov_exp(i) = 174.5;
        nacl_recovered(i) = 174.5;
    end
    if exinfo(list_nacl(i)).t_recov >= 0 && ~isempty(exinfo(list_nacl(i)).recov_fname)
        nacl_recov_exp(i) = exinfo(list_nacl(i)).id;
        if exinfo(list_nacl(i)).recov_p >= 0.05
            nacl_recovered(i) = exinfo(list_nacl(i)).id;
        end
    end
end
nacl_recov_exp(isnan(nacl_recov_exp)) = [];
nacl_recovered(isnan(nacl_recovered)) = [];
disp(['NaCl units: recovery experiment; ' num2str(length(unique(nacl_recov_exp)))]);
disp(['NaCl units: full recovery; ' num2str(length(unique(nacl_recovered)))]);

nacl_recov_exp_alone = nacl_recov_exp;
for i = 1:length(nacl_recov_exp)
    for e = 1:length(exinfo)
        if exinfo(e).id==nacl_recov_exp(i) && strcmp(exinfo(e).drugname, '5HT')
            nacl_recov_exp_alone(i) = nan;
        end
    end
end
nacl_recov_exp_alone(isnan(nacl_recov_exp_alone)) = [];
disp(['NaCl units: recovery experiments but NaCl alone; ' num2str(length(nacl_recov_exp_alone))]);

list_ser_with_recover = [];
list_nacl_recover = [];
for k = 1:length(list_ser_recover)
    for i = 1:len_ex
        if ismember(exinfo(i).id, exinfo(list_ser_recover(k)).id) && strcmp(exinfo(i).drugname, 'NaCl')
            list_nacl_recover = [list_nacl_recover, i];
            list_ser_with_recover = [list_ser_with_recover, list_ser_recover(k)];
        end
    end
end
list_nacl_recover = unique(list_nacl_recover);
list_ser_with_recover = unique(list_ser_with_recover);

    


%% visualize