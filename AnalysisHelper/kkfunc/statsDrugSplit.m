function statsDrugSplit(dat)

disp('+++++++++++++++++++++++++++++')

% 5HT
disp(['5HT: length = ' num2str(length(find(dat.is5HT==1)))])
x = dat.x(dat.is5HT==1);
y = dat.y(dat.is5HT==1);
delta1 = x - y;

[~,p] = ttest(x,y);
disp(['5HT: ttest2; p = ' num2str(p)])

p = signrank(x,y);
disp(['5HT: signrank; p = ' num2str(p)])

% NaCl
disp(['NaCl: length = ' num2str(length(find(dat.is5HT==0)))])
x = dat.x(dat.is5HT==0);
y = dat.y(dat.is5HT==0);
delta2 = x - y;

[~,p] = ttest(x,y);
disp(['NaCl: ttest2; p = ' num2str(p)])

p = signrank(x,y);
disp(['NaCl: signrank; p = ' num2str(p)])

disp('(5HT - base) - (NaCl - base)')
[~,p] = ttest2(delta1, delta2);
disp(['ttest2; p =  ' num2str(p)])

p = ranksum(delta1, delta2);
disp(['ranksum; p =  ' num2str(p)])

disp('+++++++++++++++++++++++++++++')