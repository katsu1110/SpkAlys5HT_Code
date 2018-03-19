function [para] = isolationSummary(exinfo, ver, allexp, varargin)
%% summarize isolation index and comparison to manual ranking
% INPUT: exinfo ... given by 'runExinfoAnalysis.m'
%             ver ... 1 or 2 (1 is recommended)
%
% written by Katsuhisa (18.17.17)
% ++++++++++++++++++++++++++++++++++++++++++++++++

switch nargin
    case 0
        error('Put exinfo as an input.')
    case 1
        ver = 1;
end

% % load list of index
try
    load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stim_cond_2007.mat')
catch
    load('/gpfs01/nienborg/group/Corinna/SharedCode/Katsu/incl_i_all_stim_cond_2007.mat')
end

if allexp==0    
    exinfo = exinfo(incl_i);
end

rank_b = arrayfun(@(x) x.spkqual_base, exinfo);
rank_d = arrayfun(@(x) x.spkqual_drug, exinfo);
[u,c] = uniquecount(rank_b);
disp('baseline ----------------------------------------------')
for i = 1:length(u)
    disp(['isolation rank: ' num2str(u(i)) ', counts: ' num2str(c(i))])
end
disp(' ----------------------------------------------')

disp('drug ----------------------------------------------')
[u,c] = uniquecount(rank_d);
med_base = median(rank_b(rank_b < 5));
for i = 1:length(u)
    disp(['isolation rank: ' num2str(u(i)) ', counts: ' num2str(c(i))])
end
disp(' ----------------------------------------------')
med_drug = median(rank_d(rank_d < 5));

% initialization
for i = 1:2    
    para.cond(i).out(1).val = [];
    para.cond(i).out(2).val = [];
    para.cond(i).out(3).val = [];
    para.cond(i).out(4).val = [];
    para.cond(i).rank = [];
    para.cond(i).snr = [];
    para.cond(i).Lratio = [];
    para.cond(i).iso_dist = [];
    para.cond(i).iso_score = [];
end
para.isi_all = [];
para.isi_fraction = nan(length(exinfo), 11);
thre = 1:0.1:2;

% loop 
for i = 1:length(exinfo)
    % baseline ======================
    % check if isolation-distance was computed
    
    % isolation distance
    if ~isnan(exinfo(i).isolation_distance)        
        para.cond(1).iso_dist = [para.cond(1).iso_dist, exinfo(i).isolation_distance(ver)];
        para.cond(1).out(3).val = [para.cond(1).out(3).val 0];
    else
        para.cond(1).out(3).val = [para.cond(1).out(3).val 1];
    end
    
    % manual sorting ranking
    if exinfo(i).spkqual_base==6
        para.cond(1).rank = [para.cond(1).rank, med_base];
    else
        para.cond(1).rank = [para.cond(1).rank, exinfo(i).spkqual_base];
    end
    
    % SNR
    if ~isempty(exinfo(i).snr)
        para.cond(1).snr = [para.cond(1).snr, exinfo(i).snr(ver)];
        para.cond(1).out(1).val = [para.cond(1).out(1).val 0];
    else
        para.cond(1).out(1).val = [para.cond(1).out(1).val 1];
    end
    
    % Lratio
    if ~isempty(exinfo(i).Lratio)
        para.cond(1).Lratio = [para.cond(1).Lratio, exinfo(i).Lratio(ver)];
        para.cond(1).out(2).val = [para.cond(1).out(2).val 0];
    else
        para.cond(1).out(2).val = [para.cond(1).out(2).val 1];
    end
    
    % isolation score
    if ~isempty(exinfo(i).isolation_score) && ~isnan(exinfo(i).isolation_score(1))
        para.cond(1).iso_score = [para.cond(1).iso_score, exinfo(i).isolation_score];
        para.cond(1).out(4).val = [para.cond(1).out(4).val 0];
    else
        para.cond(1).out(4).val = [para.cond(1).out(4).val 1];
    end
    
    % ISI
    if ~isempty(exinfo(i).isi)
        para.isi_all = [para.isi_all; exinfo(i).isi];
        for s = 1:11
            para.isi_fraction(i,s) = 100*length(find(exinfo(i).isi <= thre(s)))...
                /length(exinfo(i).isi);
        end
%         para.cond(1).fraction0 = [para.cond(1).fraction0; exinfo(i).isi_fraction(1,:)];
%         para.cond(1).fraction1 = [para.cond(1).fraction1; exinfo(i).isi_fraction(2,:)];
    end
    
    % drug ==========================
    % check if isolation-distance was computed
    
    % isolation distance
    if ~isnan(exinfo(i).isolation_distance_drug)        
        para.cond(2).iso_dist = [para.cond(2).iso_dist, exinfo(i).isolation_distance_drug(ver)];
        para.cond(2).out(3).val = [para.cond(2).out(3).val 0];
    else
        para.cond(2).out(3).val = [para.cond(2).out(3).val 1];
    end
    
    % manual sorting ranking
    if exinfo(i).spkqual_drug==6
        para.cond(2).rank = [para.cond(2).rank, med_drug];
    else
        para.cond(2).rank = [para.cond(2).rank, exinfo(i).spkqual_drug];
    end
    
    % SNR
    if ~isempty(exinfo(i).snr_drug)
        para.cond(2).snr = [para.cond(2).snr, exinfo(i).snr_drug(ver)];
        para.cond(2).out(1).val = [para.cond(2).out(1).val 0];
    else
        para.cond(2).out(1).val = [para.cond(2).out(1).val 1];
    end
    
    % Lratio
    if ~isempty(exinfo(i).Lratio_drug)
        para.cond(2).Lratio = [para.cond(2).Lratio, exinfo(i).Lratio_drug(ver)];
        para.cond(2).out(2).val = [para.cond(2).out(2).val 0];
    else
        para.cond(2).out(2).val = [para.cond(2).out(2).val 1];
    end
    
    % isolation score
    if ~isempty(exinfo(i).isolation_score_drug)&& ~isnan(exinfo(i).isolation_score_drug(1))
        para.cond(2).iso_score = [para.cond(2).iso_score, exinfo(i).isolation_score_drug];
        para.cond(2).out(4).val = [para.cond(2).out(4).val 0];
    else
        para.cond(2).out(4).val = [para.cond(2).out(4).val 1];
    end
    
%     % ISI fraction
%     if ~isempty(exinfo(i).isi_fraction_drug)
%         para.cond(2).fraction0 = [para.cond(2).fraction0; exinfo(i).isi_fraction_drug(1,:)];
%         para.cond(2).fraction1 = [para.cond(2).fraction1; exinfo(i).isi_fraction_drug(2,:)];
%     end
end

% % the number of units not analyzed
% leni = length(exinfo);
% disp(['baseline; no waveform stored: ' num2str(100*length(find([para.cond(1).out(1).val]==1))/leni) ' %']);
% disp(['drug; no waveform stored: ' num2str(100*length(find([para.cond(2).out(1).val]==1))/leni) ' %']);
% disp(['baseline; data with isolation distance (c0 > c1): ' ...
%     num2str(100*length(find([para.cond(1).out(3).val]==0))/length(find([para.cond(1).out(1).val]==0))) ' %']);
% disp(['drug; data with isolation distance (c0 > c1): ' ...
%     num2str(100*length(find([para.cond(2).out(3).val]==0))/length(find([para.cond(2).out(1).val]==0))) ' %']);

% visualize
close all;
h(1) = figure(1);
h(2) = figure(2);
for i = 1:4
    switch i
        case 1
            xb = para.cond(1).snr;
            yb = para.cond(1).rank(para.cond(1).out(1).val==0);
            xd = para.cond(2).snr;
            yd = para.cond(2).rank(para.cond(2).out(1).val==0);
            lab = 'S/N';
        case 2
            xb = para.cond(1).Lratio;
            yb = para.cond(1).rank(para.cond(1).out(2).val==0);
            xd = para.cond(2).Lratio;
            yd = para.cond(2).rank(para.cond(2).out(2).val==0);
            lab = 'Lratio';
        case 3
            xb = para.cond(1).iso_dist;
            yb = para.cond(1).rank(para.cond(1).out(3).val==0);
            xd = para.cond(2).iso_dist;
            yd = para.cond(2).rank(para.cond(2).out(3).val==0);
            lab = 'isolation distance';
        case 4
            xb = para.cond(1).iso_score;
            yb = para.cond(1).rank(para.cond(1).out(4).val==0);
            xd = para.cond(2).iso_score;
            yd = para.cond(2).rank(para.cond(2).out(4).val==0);
            lab = 'isolation score';
    end
    
    % histogram of spike-isolation index
    figure(1);
    subplot(4,2,2*i-1)
    histogram(xb)
    med = nanmedian(xb);
    ave = nanmean(xb);
    title(['median: ' num2str(med) ', mean: ' num2str(ave)])
    ylabel(lab)
    hold on;
    yy = get(gca,'YLim');
    plot(med*ones(1,2), yy, '-r')
    hold on;
    plot(ave*ones(1,2),yy, '-b')
%     text(med*1.05, yy(2)*0.8, ['median = ' num2str(med)], 'color', 'r')
%     text(ave*1.05, yy(2)*0.9, ['mean = ' num2str(ave)], 'color', 'b')
    set(gca, 'box','off')
    set(gca, 'TickDir','out')    
    if i == 4
        xlabel('baseline')
    end
    
    subplot(4,2,i*2)
    histogram(xd)
    med = nanmedian(xd);
    ave = nanmean(xd);    
    title(['median: ' num2str(med) ', mean: ' num2str(ave)])
    hold on;
    yy = get(gca,'YLim');
    plot(med*ones(1,2), yy, '-r')
    hold on;
    plot(ave*ones(1,2),yy, '-b')
%     text(med*1.05, yy(2)*0.8, ['median = ' num2str(med)], 'color', 'r')
%     text(ave*1.05, yy(2)*0.9, ['mean = ' num2str(ave)], 'color', 'b')
    set(gca, 'box','off')
    set(gca, 'TickDir','out')    
    if i == 4
        xlabel('drug')
    end
        
    % correlation with the manual rank
    figure(2);
    subplot(4,2,2*i -1)
    scatter(yb,xb, 20, 'filled', 'markerfacecolor','b','markerfacealpha',0.4), lsline
    [r,p] = corr(xb', yb', 'type','Spearman');
    title(['r = ' num2str(r) ', p = ' num2str(p)])
    if i == 4
        xlabel('rank (baseline)')
    end
%     xlim([0.8 4.8])
    ylabel(lab)
%     axis square
    set(gca, 'box','off')
    set(gca, 'TickDir','out')
    hold on;
    
    subplot(4,2,2*i)
    scatter(yd, xd, 20, 'filled', 'markerfacecolor','b','markerfacealpha',0.4), lsline
    [r,p] = corr(xd', yd', 'type','Spearman');
    title(['r = ' num2str(r) ', p = ' num2str(p)])
    if i == 4
        xlabel('rank (drug)')
    end    
%     xlim([0.8 4.8])
%     axis square    
    set(gca, 'box','off')
    set(gca, 'TickDir','out')    
    hold on;
end
set(h(1), 'position', [543   140   560   800])
set(h(2), 'position', [1120         141         517         799])

% ISI fraction visualization
h(3) = figure(3);
subplot(3,4,1)
histogram(para.isi_all)
title('ISI of c1 (ms)')
yy = get(gca, 'YLim');
med = nanmedian(para.isi_all);
ave = nanmean(para.isi_all);
hold on;
plot(med*ones(1,2), yy, '-r')
text(med+0.25, yy(2)*0.8, ['median = ' num2str(med)], 'color', 'r')
hold on;
plot(ave*ones(1,2), yy, '-b')
text(ave+0.25, yy(2)*0.9, ['mean = ' num2str(ave)], 'color', 'b')
set(gca, 'box', 'off')
set(gca, 'TickDir', 'out')
axis square
for i = 1:11
    subplot(3,4,i+1)
    histogram(para.isi_fraction(:,i));
    title(['cutoff = ' num2str(thre(i)) ' ms'])
    yy = get(gca, 'YLim');
    less2 = sum(para.isi_fraction(:,i) < 2);
    more2 = sum(para.isi_fraction(:,i) >= 2);
    med = nanmedian(para.isi_fraction(:,i));
    ave = nanmean(para.isi_fraction(:,i));
    hold on;
    plot(med*ones(1,2), yy, '-r')
    text(med+0.05, yy(2)*0.8, ['median = ' num2str(med)], 'color', 'r')
    hold on;
    plot(ave*ones(1,2), yy, '-b')
    text(ave+0.05, yy(2)*0.9, ['mean = ' num2str(ave)], 'color', 'b')   
    text(2, yy(1) + (yy(2)-yy(1))*0.25, ['< 2%: ' num2str(less2)],'color','m')
    text(2, yy(1) + (yy(2)-yy(1))*0.45, ['>= 2%: ' num2str(more2)],'color','g')
    set(gca, 'box', 'off')
    set(gca, 'TickDir', 'out')
    axis square
    if i==5
%         ylabel('ISI fraction: c1 [%]')
        ylabel('counts')
    end
    if i==9
        xlabel('ISI fraction: c1 [%]')
    end
end
set(gcf, 'position', [680   420   786   558])

% compare included and excluded experiments
if allexp==1
    h(4) = figure(4);
    v = 1:length(exinfo);
    for i = 1:4
        vb = v;
        vd = v;
        switch i
            case 1
                xb = para.cond(1).snr;
                vb(para.cond(1).out(1).val==1) = [];
                xb_inc = log(xb(ismember(vb, incl_i)));
                xb_exc = log( xb(~ismember(vb, incl_i)));
                xd = para.cond(2).snr;
                vd(para.cond(2).out(1).val==1) = [];
                xd_inc = log(xd(ismember(vd, incl_i)));
                xd_exc = log(xd(~ismember(vd, incl_i)));
                xrange = [0.01 4];
                lab = 'S/N (log)';
            case 2
                xb = para.cond(1).Lratio;
                vb(para.cond(1).out(2).val==1) = [];
                xb_inc = log(xb(ismember(vb, incl_i)));
                xb_exc = log(xb(~ismember(vb, incl_i)));
                xd = para.cond(2).Lratio;
                vd(para.cond(2).out(2).val==1) = [];
                xd_inc = log(xd(ismember(vd, incl_i)));
                xd_exc = log(xd(~ismember(vd, incl_i)));
                xrange = [-50 -0.01];
                lab = 'Lratio (log)';
            case 3
                xb = para.cond(1).iso_dist;
                vb(para.cond(1).out(3).val==1) = [];
                xb_inc = log(xb(ismember(vb, incl_i)));
                xb_exc = log(xb(~ismember(vb, incl_i)));
                xd = para.cond(2).iso_dist;
                vd(para.cond(2).out(3).val==1) = [];
                xd_inc = log(xd(ismember(vd, incl_i)));
                xd_exc = log(xd(~ismember(vd, incl_i)));
                xrange = [0.01 15];
                lab = 'isolation distance (log)';
            case 4
                xb = para.cond(1).iso_score;
                vb(para.cond(1).out(4).val==1) = [];
                xb_inc = xb(ismember(vb, incl_i));
                xb_exc = xb(~ismember(vb, incl_i));
                xd = para.cond(2).iso_score;
                vd(para.cond(2).out(4).val==1) = [];
                xd_inc = xd(ismember(vd, incl_i));
                xd_exc = xd(~ismember(vd, incl_i));
                xrange = [0.3 1];
                lab = 'isolation score';
        end
        
        disp([lab '; excluded:  ' num2str(length(xb_exc))])
        disp([lab '; included:  ' num2str(length(xb_inc))])
        
        subplot(2,4,i)
        histogram(xb_exc)
        hold on;
        histogram(xb_inc)
        xlim(xrange)
        title(lab)
        set(gca, 'box', 'off')
        set(gca, 'TickDir', 'out')
        if i==1
            ylabel('baseline')
        end
        
        subplot(2,4,i+4)
        histogram(xd_exc)
        hold on;
        histogram(xd_inc)
        xlim(xrange)
        set(gca, 'box', 'off')
        set(gca, 'TickDir', 'out')
        if i==1
            ylabel('drug')
        end
        
    end
    set(gcf, 'position', [ 680   558   876   420])
end