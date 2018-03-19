function b = plotPuGLM(exinfo, drugtype, stmtype)

% specify stmtype
if strcmp(stmtype, 'rc')
    load('Z:\Corinna\SharedCode\Katsu\list_RC.mat')
    exinfo = exinfo(list_RC);
    if drugtype==1
        exinfo = exinfo([exinfo.is5HT]==1);
    else
        exinfo = exinfo([exinfo.is5HT]==0);
    end
elseif strcmp(stmtype, 'all')
%     load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stimulus_cond.mat')
%     exinfo = exinfo(incl_i);
    load('Z:\Corinna\SharedCode\Katsu\list_unit.mat')
    exinfo = exinfo(list_unit);    
    if drugtype==1
        exinfo = exinfo([exinfo.is5HT]==1);
    else
        exinfo = exinfo([exinfo.is5HT]==0);
    end
    
else
    load('Z:\Corinna\SharedCode\Katsu\list_unit.mat')
    exinfo = exinfo(list_unit);
    exinfo = exinfo([exinfo.is5HT]==1);
    len_exp = length(exinfo);
    match = zeros(1, len_exp);
    for i = 1:len_exp
        if strcmp(exinfo(i).param1, stmtype)
            match(i) = 1;
        end
    end
    exinfo(match==0) = []; 
end

len_exp = length(exinfo);
disp(['the number of experiments for ' stmtype ': ' num2str(len_exp)])

bdrug = nan(len_exp,1);
bps = nan(len_exp,1);
bintr = nan(len_exp,1);
% bdrug_add = nan(len_exp,1);
% bps_add = nan(len_exp,1);
% bintr_add = nan(len_exp,1);

close all;
figure;


% pupil size time-course
subplot(2,2,1)
ncol = 1000;
for i = 1:len_exp
    s = size(exinfo(i).psglm.tc,2);
    if ncol > s
        ncol = s;
    end
end
mat = [];
for i = 1:len_exp
    mat = [mat; exinfo(i).psglm.tc(:,1:ncol)];
end
nstm = unique(mat(:,2));
nlen = length(nstm);
if nlen == 1
    offset = 0;
else 
    offset = 0.1;
end
ncol = ncol - 2;
for d = 1:2
    switch d
        case 1
            col = [0 0 0];
        case 2
            col = [1 0 0];
    end
    for n = 1:length(nstm)
        temp = mat(mat(:,2)==nstm(n) & mat(:,1)==d-1, 3:end);
        me = mean(temp, 1);
        sem = std(temp, [], 1)/sqrt(size(temp,1));
        fill_between(([[1:ncol] + ncol*(n-1)]/500) + offset*(n-1),...
            me - sem, me + sem, col)
        hold on;
        plot(([[1:ncol] + ncol*(n-1)]/500) + offset*(n-1), me, '-','color',col)
        hold on;
    end
end
xlabel('time (s)')
ylabel('pupil size (a.u.)')
set(gca,'box','off')
set(gca,'TickDir','out')
axis square

% interaction table
subplot(2,2,2)
intrmat = zeros(len_exp,4);
for i = 1:len_exp
    intrmat(i,:) = reshape(exinfo(i).psglm.inter_table, 1, 4);
end
intr_me = reshape(mean(intrmat,1),2,2);
% intr_sem = reshape(std(intrmat, [], 1)/sqrt(len_exp), 2,2);
imagesc(1:2, 1:2, intr_me);
colormap(jet);
colorbar('southoutside')
hold on;
plot([0.5 2.5],[1.5 1.5],'-w','linewidth',2)
hold on;
plot([1.5 1.5],[0.5 2.5],'-w','linewidth',2)
hold on;
text(0.7, 1, [num2str(intr_me(1,1)-1)],'color','w')
text(0.7, 2, [num2str(intr_me(2,1)-1)],'color','w')
text(1.7, 1, [num2str(intr_me(1,2)-1)],'color','w')
text(1.7, 2, [num2str(intr_me(2,2)-1)],'color','w')
% text(0.7, 1, [num2str(intr_me(2,1)-1) ' ± ' num2str(intr_sem(2,1))])
% text(0.7, 2, [num2str(intr_me(1,1)-1) ' ± ' num2str(intr_sem(1,1))])
% text(1.7, 1, [num2str(intr_me(2,2)-1) ' ± ' num2str(intr_sem(2,2))])
% text(1.7, 2, [num2str(intr_me(1,2)-1) ' ± ' num2str(intr_sem(1,2))])
axis([0.5 2.5 0.5 2.5])
xlabel('drug condition')
ylabel('pupil size')
set(gca,'XTick',[1 2],'XTickLabel',{'5HT','baseline'})
set(gca,'YTick',[1 2],'YTickLabel',{'small','large'})
set(gca,'box','off')
set(gca,'TickDir','out')
axis square
xtickangle(45)


% model performance
subplot(2,2,3)
plot([0.5 4.5],[0 0],'-k')
perf = nan(len_exp, 4);
out = zeros(1, len_exp);
for i = 1:len_exp
    if min(exinfo(i).psglm.perf) > 0 && max(exinfo(i).psglm.perf) < 1
        perf(i,:) = exinfo(i).psglm.perf(1:4);
    else
        out(i) = 1;
    end
end
waterfallchart(nanmean(perf, 1));
hold on;
errorbar(1:4, nanmean(perf, 1), nanstd(perf, [], 1)/sqrt(sum(out==0)),'.k')
xlim([0.5 4.5])
xlabel('predictors')
ylabel('variance explained')
for i = 1:3
    [~,p] = ttest(perf(:,i), perf(:,i+1));
    text(0.6 + i, 0.8, ['p = ' num2str(p)])
end
ylim([min(nanmean(perf,1))-0.05 max(nanmean(perf,1))+0.1]);
set(gca,'XTick',1:4,'XTickLabel',{'stm','+drug','+pupil','+drug x pupil'})
set(gca,'box','off')
set(gca,'TickDir','out')
xtickangle(45)

% beta weights
subplot(2,2,4)
plot([0.5 3.5],[0 0],'-k')
hold on;
for i = 1:len_exp
    if out(i)==0
        bdrug(i) = exinfo(i).psglm.b(3,1);
        bps(i) = exinfo(i).psglm.b(3,2);
        bintr(i) = exinfo(i).psglm.b(3,3);
%         bdrug_add(i) = exinfo(i).psglm.b(6,4);
%         bps_add(i) = exinfo(i).psglm.b(6,5);
%         bintr_add(i) = exinfo(i).psglm.b(6,6);

        scatter(1:3, [bdrug(i) bps(i) bintr(i)], 30, 'filled', 'markerfacecolor', 'r', ...
        'markeredgecolor','r','markerfacealpha',0.4,'markeredgealpha',0.8)
        hold on;
    end
end
errorbar(1:3,[nanmean(bdrug) nanmean(bps) nanmean(bintr)],...
    [nanstd(bdrug) nanstd(bps) nanstd(bintr)]/sqrt(sum(out==0)),'color',[0 0.5 0])

xlim([0.5 3.5])
set(gca,'XTick',1:3,'XTickLabel',{'drug','pupil','drug x pupil'})
set(gca,'box','off')
set(gca,'TickDir','out')
xtickangle(45)
xlabel('predictors')
ylabel('beta weights (a.u.)')

b = [bdrug, bps, bintr];
for i = 1:3
    p = signrank(b(:,i));
    text(i - 0.15, 0.1, ['p = ' num2str(p)])
end

set(gcf, 'position', [346   367   894   278])

