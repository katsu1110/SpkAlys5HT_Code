function [unit, list_new, animal] = list2unit(exinfo, list)
%% convert list to unit id and output the list of rows including one unique unit
% list_new contains list of rows without overlap of units (only the first
% experiment was included)

unit = [exinfo(list).id];

% unique count
u = unique(unit);
lenu = length(u);
c = ones(1, lenu);
for k = 1:lenu
    c(k) = sum(unit==u(k));
end

% remove overlapped units
list_new = nan(1,lenu);
animal = zeros(1, lenu);
for k = 1:lenu
    row = find(unit==u(k));
    animal(k) = exinfo(list(row(1))).ismango;
    if c(k)==1
        list_new(k) = list(row);
    else
        disp('+++++++++++++++++++++++')
        disp(['overlapped rows: ' num2str(list(row))])
        disp('+++++++++++++++++++++++')
        min_t = 10000;
        remain = 1;
        for r = 1:c(k)
            f = exinfo(list(row(r))).fname_drug;
            disp(['unit ' num2str(u(k)) ': ' f])
            under = strfind(f,'_');
            grating = strfind(f,'grating');
            cand = f(under(end)+1:grating(1)-2);
            p = strfind(cand,'PM');
            if ~isempty(p)
                cand = cand(1:p-1);
                cand = [num2str(12 + str2double(cand(1))) cand(2:end)];
            end
            comma = strfind(cand,'.');
            time = 60*str2double(cand(1:comma-1)) + str2double(cand(comma+1:end));
            if min_t > time
                min_t = time;
                remain = r;
            end
        end
        disp('-----------------------------')
        disp(['earlyest: ' exinfo(list(row(remain))).fname_drug])
        disp('-----------------------------')
        list_new(k) = list(row(remain));
    end
end
               