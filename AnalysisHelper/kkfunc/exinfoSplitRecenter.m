function [exinfo0, exinfo1] = exinfoSplitRecenter(exinfo)
% split exinfo based on the number of recenterings

recentering = nan(1, length(exinfo));
for i = 1:length(exinfo)
    val_base = size(exinfo(i).recentering_win,1) -1;
    val_drug = size(exinfo(i).recentering_win_drug, 1) - 1;
    if val_base==0 && val_drug==0
      recentering(i) = 0;
    else
      recentering(i) = 1;
    end
end

exinfo0 = exinfo(recentering==0);
exinfo1 = exinfo(recentering==1);