function plotIsoDistvsFeatures(exinfo, drug, isolationIndex)

% lista = [3     4    23    35    36 ...
%     38    44    57    85    88 ...
%     91   104   105   106   124 ...
%    145   152];

close all
len_ex = length(exinfo);

switch isolationIndex
    case 1
        name = 'isolation_distance';
    case 2
        name = 'Lratio';
end

features = {'energy + PC1', '+ PC2', '+ PC3', '+ time', '+ peak-valley',...
             '+ peak', '+ valley', '+ area', '+ non-linear energy', '+ PC4', '+ PC5', ...
             '+ PC6', '+ PC7', '+ PC8'};
         
N = length(features);

for l = 1:len_ex

    isodist = nan(1, N);
    featc = 1;
    for i = 1:N
        featc = [featc i+1];
        [spkiso] = spikeisolation(exinfo(l), drug, featc, 0,0);
        isodist(i) = spkiso.(name)(2);
    end

    if drug==0
        f = exinfo(l).fname;
    else
        f = exinfo(l).fname_drug;
    end

    % close all;
    figure;
    bar(1:N, isodist)
    ylabel('isolation distance')
    title(f)
    set(gca, 'XTick', 1:N, 'XTickLabel', features)
    xtickangle(45)
    set(gca, 'box', 'off')
    set(gca, 'TickDir', 'out')

end