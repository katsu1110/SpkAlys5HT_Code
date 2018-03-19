function [spkiso, h] = spikeisolation(exinfo, drug, featc, fig, saveoption, varargin)
%% compute spike-isolation quality (S/N, ISI, L, Lratio, islation distance, isolation score)
% INPUT: fname ... directory of ex-file with c1 or c2 unit with recording data
%             drug ... 0, baseline 1, drug condition
%             featc ... vector specifying columns of features to use (energy, PC1 - PC3, time, peak - valley,
%             peak, valley, area, non-linear energy, PC4 - PC8) 
%             exp, [2 3 5]
%             fig ... if 1, plot results
%             saveoption ... 1 save fig
%
% written by Katsuhisa (14.07.17)
% +++++++++++++++++++++++++++++++++++++++++++++++++

switch nargin
    case 0
        error('ex-file with sorted recording data must be given.')
    case 1
        drug = 0;
        featc = [1 2];
        fig = 0;
        saveoption = 0;
    case 2
        featc = [1 2];
        fig = 0;
        saveoption = 0;
    case 3
        fig = 0;
        saveoption = 0;
    case 4
        saveoption = 0;
end

if featc > 14
    featc = 14;
elseif featc < 1
    disp('the number of features is adjusted to 2.')
    featc = [1 2];
end

if drug==0
    fname = exinfo.fname;
    rank = exinfo.spkqual_base;
    postfix = '_spkiso_base';
else
    fname = exinfo.fname_drug;
    rank = exinfo.spkqual_drug;
    postfix = '_spkiso_drug';
end

% load c1 or c2
load(fname)
ex1 = ex;

% load c0
under = strfind(fname, '_');
fname(under(2):under(2)+2) = '_c0';
load(fname)
ex0 = ex;

tr  = find(abs([ex.Trials.Reward])>0);
    
% struct initialization
spkiso = struct('c0',[],'c1',[],'snr',[nan nan],'fraction', nan(2,11), 'L',[nan nan],'Lratio',[nan nan],'isolation_distance',[nan nan],'isolation_score',nan);

% check the length of spike wave and take the least, stimlus duration
l = 1000;
t = nan(1, length(tr));
for i = 1:length(tr)
    % shape of waveform
    wl0 = size(ex0.Trials(tr(i)).Waves,2);
    wl1 = size(ex1.Trials(tr(i)).Waves,2);
    
    wl = min([wl0, wl1]);
    if wl < l
        l = wl;
    end
    
    % stimulus duration
    t(i) = max(ex.Trials(tr(i)).Start - ex.Trials(tr(i)).TrialStart);
end
t = round(100*median(t))/100;

% compute spike-isolation index for each spike
spkiso.c0.time = [];
spkiso.c1.time = [];
spkiso.c0.refp = [];
spkiso.c1.refp = [];
spkiso.c0.area = [];
spkiso.c1.area = [];
% spkiso.c0.fwhm = [];
% spkiso.c1.fwhm = [];
spkiso.c0.energy = [];
spkiso.c1.energy = [];
spkiso.c0.nlenergy = [];
spkiso.c1.nlenergy = [];
spkiso.c0.waves = [];
spkiso.c1.waves = [];
spkiso.c0.norm_waves = [];
spkiso.c1.norm_waves = [];
for i = 1:length(tr)
    % get spikes within the stimulus presentation
    c0spk = ex0.Trials(tr(i)).Spikes;
    c1spk = ex1.Trials(tr(i)).Spikes;
    spkt0 = find(c0spk > 0 & c0spk <= t);
    spkt1 = find(c1spk > 0 & c1spk <= t);
    c0spk = c0spk(spkt0);
    c1spk = c1spk(spkt1);
    
    c0wave = ex0.Trials(tr(i)).Waves;
    c1wave = ex1.Trials(tr(i)).Waves;
    c0wave = c0wave(spkt0, :);
    c1wave = c1wave(spkt1, :);      
            
    % time
    stime = ex.Trials(tr(i)).TrialStart - ex.Trials(1).TrialStart;
    spkiso.c0.time = [spkiso.c0.time; c0spk + stime];
    spkiso.c1.time = [spkiso.c1.time; c1spk + stime];
    
    % refractory period
    spkiso.c0.refp = [spkiso.c0.refp; diff(c0spk)]; 
    spkiso.c1.refp = [spkiso.c1.refp; diff(c1spk)]; 
    
    % energy, store amplitude of waves, store waves    
    for k = 1:length(spkt0)
        [feature] = featureCalc(c0wave(k, 1:l));
        spkiso.c0.area = [spkiso.c0.area; feature.area];
%         spkiso.c0.fwhm = [spkiso.c0.fwhm; feature.fwhm];
        spkiso.c0.energy = [spkiso.c0.energy; feature.energy];
        spkiso.c0.nlenergy = [spkiso.c0.nlenergy; feature.nlenergy];
        spkiso.c0.waves = [spkiso.c0.waves; c0wave(k,1:l)];
        spkiso.c0.norm_waves = [spkiso.c0.norm_waves; c0wave(k,1:l)/feature.energy];
    end
    for k = 1:length(spkt1)
        [feature] = featureCalc(c1wave(k, 1:l));
        spkiso.c1.area = [spkiso.c1.area; feature.area];
%         spkiso.c1.fwhm = [spkiso.c1.fwhm; feature.fwhm];
        spkiso.c1.energy = [spkiso.c1.energy; feature.energy];
        spkiso.c1.nlenergy = [spkiso.c1.nlenergy; feature.nlenergy];
        spkiso.c1.waves = [spkiso.c1.waves; c1wave(k,1:l)];
        spkiso.c1.norm_waves = [spkiso.c1.norm_waves; c1wave(k,1:l)/feature.energy];
    end
end

% other features
% spkiso.c0.peak = max(spkiso.c0.waves, [], 2);
% spkiso.c0.valley = min(spkiso.c0.waves, [], 2);
% spkiso.c1.peak = max(spkiso.c1.waves, [], 2);
% spkiso.c1.valley = min(spkiso.c1.waves, [], 2);

spkiso.c0.peak = max(spkiso.c0.norm_waves, [], 2);
spkiso.c0.valley = min(spkiso.c0.norm_waves, [], 2);
spkiso.c1.peak = max(spkiso.c1.norm_waves, [], 2);
spkiso.c1.valley = min(spkiso.c1.norm_waves, [], 2);

% S/N ----------------------------------------------------------
% meanwave0 = mean(spkiso.c0.waves, 1);
meanwave1 = mean(spkiso.c1.waves, 1);
resid1 = spkiso.c1.waves - meanwave1;
C = 5;
spkiso.snr(1) = (max(meanwave1) - min(meanwave1))/C*std(resid1(:));
spkiso.snr(2) = var(spkiso.c1.waves(:))/var(spkiso.c0.waves(:));


% ISI ---------------------------------------------------------
thre = 1:0.1:2;
for i = 1:11
    % c0
    spkiso.fraction(1,i) = 100*length(find(spkiso.c0.refp < thre(i)*0.001))/...
        length(spkiso.c0.refp);
    
    % c1
    spkiso.fraction(2,i) = 100*length(find(spkiso.c1.refp < thre(i)*0.001))/...
        length(spkiso.c1.refp);
end


% PCA for the spike waveform --------------------
len_spk0 = size(spkiso.c0.waves,1);
len_spk1 = size(spkiso.c1.waves,1);

X = [spkiso.c0.waves; spkiso.c1.waves];
sigma = (X'*X)/size(X,1);
[U,S,V] = svd(sigma);
pca = X*U(:,1:8);

% c0
spkiso.c0.pca = pca(1:len_spk0,:);

% c1
spkiso.c1.pca = pca(len_spk0+1:end,:);

X = [spkiso.c0.norm_waves; spkiso.c1.norm_waves];
sigma = (X'*X)/size(X,1);
[U,S,V] = svd(sigma);
pca = X*U(:,1:8);

% c0
spkiso.c0.norm_pca = pca(1:len_spk0,:);

% c1
spkiso.c1.norm_pca = pca(len_spk0+1:end,:);


% The Mahalanobis distance ------------------------
% c0
X0 = [spkiso.c0.energy, spkiso.c0.norm_pca(:,1)];
X00 = [spkiso.c0.energy, spkiso.c0.norm_pca(:,1:3), spkiso.c0.time, ...
    spkiso.c0.peak - spkiso.c0.valley, spkiso.c0.peak, spkiso.c0.valley, ...
    spkiso.c0.area, spkiso.c0.nlenergy, spkiso.c0.norm_pca(:,4:8)];

% c1
X1 = [spkiso.c1.energy, spkiso.c1.norm_pca(:,1)];
X11 = [spkiso.c1.energy, spkiso.c1.norm_pca(:,1:3), spkiso.c1.time, ...
    spkiso.c1.peak - spkiso.c1.valley, spkiso.c1.peak, spkiso.c1.valley, ...
    spkiso.c1.area, spkiso.c1.nlenergy, spkiso.c1.norm_pca(:,4:8)];

% feature selection
X00 = X00(:, featc);
X11 = X11(:, featc);

% % zscoring
% X00 = zscore(X00, 0, 1);
% X11 = zscore(X11, 0, 1);

% mahalanobis for c0
% within
spkiso.c0.mahalanobis(:,1) = mahal(X0, X0);
spkiso.c0.mahalanobis2(:,1) = mahal(X00, X00);

% outside
spkiso.c0.mahalanobis(:,2) = mahal(X0, X1);
spkiso.c0.mahalanobis2(:,2) = mahal(X00, X11);

% mahalanobis for c1
% within
spkiso.c1.mahalanobis(:,1) = mahal(X1, X1);
spkiso.c1.mahalanobis2(:,1) = mahal(X11, X11);

% outside
spkiso.c1.mahalanobis(:,2) = mahal(X1, X0);
spkiso.c1.mahalanobis2(:,2) = mahal(X11, X00);


% Lratio -----------------------------------------------------------
spkiso.L = [0 0];
len_feat = length(featc);
for i = 1:len_spk0
    spkiso.L(1) = spkiso.L(1) + (1 - chi2cdf(spkiso.c0.mahalanobis(i,2),2));
    spkiso.L(2) = spkiso.L(2) + (1 - chi2cdf(spkiso.c0.mahalanobis2(i,2), len_feat));
end
spkiso.Lratio = spkiso.L/len_spk1;

% Isolation distance ----------------------------------------
if len_spk0 > len_spk1
    sorted_md = sort(spkiso.c0.mahalanobis(:,2));
    spkiso.isolation_distance(1) = sorted_md(len_spk1);
    sorted_md = sort(spkiso.c0.mahalanobis2(:,2));
    spkiso.isolation_distance(2) = sorted_md(len_spk1);
else
    iso_dist = nan(2,10);
    seed = [1220, 0217, 1098454, 34392, 223421, 28942, ...
        394859, 12837463, 22385, 985738, 2736482];
    for r = 1:10
        
        % random sampling
        rng(seed(r))
        idx = datasample(1:len_spk1, len_spk0, 'Replace', false);
        
        % PCA for the spike waveform --------------------
        X = [spkiso.c0.norm_waves; spkiso.c1.norm_waves(idx,:)];
        sigma = (X'*X)/size(X,1);
        [U,S,V] = svd(sigma);
        pca = X*U(:,1:8);

        % c0
        norm_pca0 = pca(1:len_spk0,:);

        % c1
        norm_pca1 = pca(len_spk0+1:end,:);

        % The Mahalanobis distance ------------------------
        % c0
        X0 = [spkiso.c0.energy, norm_pca0(:,1)];
        X00 = [spkiso.c0.energy, spkiso.c0.norm_pca(:,1:3), spkiso.c0.time, ...
            spkiso.c0.peak - spkiso.c0.valley, spkiso.c0.peak, spkiso.c0.valley, ...
            spkiso.c0.area, spkiso.c0.nlenergy, spkiso.c0.norm_pca(:,4:8)];

        % c1
        X1 = [spkiso.c1.energy(idx), norm_pca1(:,1)];
        X11 = [spkiso.c1.energy(idx), spkiso.c1.norm_pca(idx,1:3), spkiso.c1.time(idx), ...
            spkiso.c1.peak(idx) - spkiso.c1.valley(idx), spkiso.c1.peak(idx), spkiso.c1.valley(idx), ...
            spkiso.c1.area(idx), spkiso.c1.nlenergy(idx), spkiso.c1.norm_pca(idx,4:8)];

        % feature selection
        X00 = X00(:, featc);
        X11 = X11(:, featc);

        % mahalanobis for c0
        % outside
        mahalanobis = mahal(X0, X1);
        mahalanobis2 = mahal(X00, X11);
        
        % isolation distance -------------------------------------
        iso_dist(1,r) = max(mahalanobis);
        iso_dist(2,r) = max(mahalanobis2);
    end
    
    spkiso.isolation_distance(1) = mean(iso_dist(1,:),2);
    spkiso.isolation_distance(2) = mean(iso_dist(2,:),2);
end

% % isolation score -----------------------------------------
% % normalization constant
% d = 0;
% for i = 1:len_spk1
%     wmat = spkiso.c1.waves;
%     wmat(i,:) = [];
%     d = d + sum(sum((wmat - spkiso.c1.waves(i,:)).^2, 2));
% end
% d0 = d/(len_spk1*(len_spk1 - 1));
% 
% % similarity index
% si = 0;
% for i = 1:len_spk1
%     wmat = spkiso.c1.waves;
%     wmat(i,:) = [];
%     sig = sum((wmat - spkiso.c1.waves(i,:)).^2, 2);
%     dd = sum(exp(-(10/d0)*[sig; ...
%         sum((spkiso.c0.waves - spkiso.c1.waves(i,:)).^2, 2)]));
%     si = si + sum(exp(-(10/d0)*sig)/dd);
% end
% spkiso.isolation_score = si/len_spk1;

if fig==1
    col0 = colormap(winter(7));
    close(gcf);
    col1 = colormap(autumn(7));
    close(gcf);
    
    bin = 5;
    c0perbin = floor(len_spk0/bin);
    c1perbin = floor(len_spk1/bin);
    
    h = figure;
    begin0 = 1;
    begin1 = 1;
    for b = 1:bin
        % PCA    
        subplot(2,3,1)
        scatter3(spkiso.c0.time(begin0:begin0+c0perbin-1), spkiso.c0.pca(begin0:begin0+c0perbin-1,1), spkiso.c0.pca(begin0:begin0+c0perbin-1,2), ...
            10, 'o', ...
            'markerfacecolor', col0(b,:), 'markeredgecolor',col0(b,:),...
            'markerfacealpha', 0.2, 'markeredgealpha', 0.8);
        hold on;
        scatter3(spkiso.c1.time(begin1:begin1+c1perbin-1), spkiso.c1.pca(begin1:begin1+c1perbin-1,1), spkiso.c1.pca(begin1:begin1+c1perbin-1,2),...
            10, 'o',...
            'markerfacecolor', col1(b,:), 'markeredgecolor',col1(b,:),...
            'markerfacealpha', 0.2, 'markeredgealpha', 0.8);
        ylabel('PC1')
        zlabel('PC2')
        set(gca,'TickDir','out')
        set(gca,'box','off')
        axis square

        % PC1 vs Energy
        subplot(2,3,2)
        p0 = scatter3(spkiso.c0.time(begin0:begin0+c0perbin-1), spkiso.c0.pca(begin0:begin0+c0perbin-1,1), spkiso.c0.energy(begin0:begin0+c0perbin-1), ...
            10, 'o', ...
            'markerfacecolor', col0(b,:), 'markeredgecolor',col0(b,:),...
            'markerfacealpha', 0.2, 'markeredgealpha', 0.8);
        hold on;
        p1 = scatter3(spkiso.c1.time(begin1:begin1+c1perbin-1), spkiso.c1.pca(begin1:begin1+c1perbin-1,1), spkiso.c1.energy(begin1:begin1+c1perbin-1), ...
             10, 'o',...
            'markerfacecolor', col1(b,:), 'markeredgecolor',col1(b,:),...
            'markerfacealpha', 0.2, 'markeredgealpha', 0.8);
        legend([p0,p1],'c0','c1')
        legend('boxoff')
        ylabel('PC1')
        zlabel('energy')
        set(gca,'TickDir','out')
        set(gca,'box','off')
        axis square

        % PCA (energy-normalized)
        subplot(2,3,4)
        scatter3(spkiso.c0.time(begin0:begin0+c0perbin-1), spkiso.c0.norm_pca(begin0:begin0+c0perbin-1, 1), spkiso.c0.norm_pca(begin0:begin0+c0perbin-1,2), ...
            10, 'o',...
            'markerfacecolor', col0(b,:), 'markeredgecolor',col0(b,:),...
            'markerfacealpha', 0.2, 'markeredgealpha', 0.8);
        hold on;
        scatter3(spkiso.c1.time(begin1:begin1+c1perbin-1),spkiso.c1.norm_pca(begin1:begin1+c1perbin-1,1), spkiso.c1.norm_pca(begin1:begin1+c1perbin-1,2), ...
             10, 'o',...
            'markerfacecolor', col1(b,:), 'markeredgecolor',col1(b,:),...
            'markerfacealpha', 0.2, 'markeredgealpha', 0.8);
        ylabel('PC1 (energy normalized)')
        zlabel('PC2 (energy normalized)')
        set(gca,'TickDir','out')
        set(gca,'box','off')
        axis square

        % PC1 (energy-normalized) vs Energy
        subplot(2,3,5)
        p0 = scatter3(spkiso.c0.time(begin0:begin0+c0perbin-1),spkiso.c0.norm_pca(begin0:begin0+c0perbin-1,1), spkiso.c0.energy(begin0:begin0+c0perbin-1), ...
             10, 'o',...
            'markerfacecolor', col0(b,:), 'markeredgecolor',col0(b,:),...
            'markerfacealpha', 0.2, 'markeredgealpha', 0.8);
        hold on;
        p1 = scatter3(spkiso.c1.time(begin1:begin1+c1perbin-1),spkiso.c1.norm_pca(begin1:begin1+c1perbin-1,1), spkiso.c1.energy(begin1:begin1+c1perbin-1), ...
         10, 'o',...
             'markerfacecolor', col1(b,:), 'markeredgecolor',col1(b,:),...
            'markerfacealpha', 0.2, 'markeredgealpha', 0.8);
        legend([p0,p1],'c0','c1')
        legend('boxoff')
        ylabel('PC1 (energy normalized)')
        zlabel('energy')
        set(gca,'TickDir','out')
        set(gca,'box','off')
        axis square
        
        begin0 = begin0 + c0perbin;
        begin1 = begin1 + c1perbin;
    end    
    
    % refractory period (inter-spike intervals)
    subplot(2,3,3)
    c0isi = spkiso.c0.refp*1000;
    c1isi = spkiso.c1.refp*1000;
    histogram(c0isi(c0isi < 5),'binwidth',0.1)
    hold on;
    histogram(c1isi(c1isi < 5),'binwidth',0.1)
    hold on;
    yy = get(gca, 'YLim');
    plot(2*ones(1,2), yy, '-r')
    hold on;
    for k = 1:11
        text(5.05, yy(1)+k*(yy(2)-yy(1))/12, ['<' num2str(thre(k)) '(ms): ' num2str(spkiso.fraction(2,k)) '%'])
    end
    xlim([0 5])
    xlabel('inter-spike intervals (ms)')
    ylabel('counts')
    set(gca,'TickDir','out')
    set(gca,'box','off')
    axis square
    
    title(['manual rank: ' num2str(rank)])

    % The Mahalanobis distance
    subplot(2,3,6)
    p0 = histogram(spkiso.c0.mahalanobis2(:,2));
    hold on;
    p1 = histogram(spkiso.c1.mahalanobis2(:,1));
    xx = get(gca, 'XLim');
    xlim([0 xx(2)])
    legend([p0,p1],'c0','c1')
    legend('boxoff')
    xlabel('The mahalanobis distance')
    ylabel('counts')
    title(['Lratio = ' num2str(spkiso.Lratio(2)) ', iso-dist = ' num2str(spkiso.isolation_distance(2))])
    set(gca,'TickDir','out')
    set(gca,'box','off')
%     set(gca, 'xscale','log')
    axis square
    
    set(gcf, 'position', [680   558   842   420])
    
    savedir = 'Z:\Katsuhisa\interaction_project\dataset\Figures\Spike_isolation\';
    if saveoption==1
        set(h,'name',fname,'NumberTitle','off')
        savefig(h, [savedir exinfo.figname postfix '.fig'])
    end
end

% subfunction to calculate area, energy, non-linear energy
function [feature] = featureCalc(wave)
% calculate features based on the spike waveform

lenw = length(wave);

% % peak 
% peak = max(wave, [], 2);
% 
% % valley
% valley = min(wave, [], 2);
% 
% % amplitude
% amp = peak - valley;
% 
% % fwhm
% width = fwhm(1:lenw, wave);

% area
area = sum(abs(wave));

% energy
energy =  sqrt(sum(wave.^2))/lenw;

% non-linear energy
temp1 = [0 wave(1:end-1)];
temp2 = [wave(2:end) 0];
nl_energy = sqrt(sum(wave.^2 - temp1.*temp2))/lenw;

% structurize
feature = struct('area', area, 'energy', energy, 'nlenergy', nl_energy);




