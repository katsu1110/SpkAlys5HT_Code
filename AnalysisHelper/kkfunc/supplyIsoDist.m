function [exinfo] = supplyIsoDist(exinfo)

for i = 1:length(exinfo)
%     if isempty(exinfo(i).isolation_distance) || isempty(exinfo(i).isolation_distance_drug) ...
%             || isnan(exinfo(i).isolation_distance(1)) || isnan(exinfo(i).isolation_distance_drug(1))
        try
            [spkiso] = isodist(exinfo(i), 0);
            exinfo(i).isolation_distance = spkiso.isolation_distance;
            [spkiso] = isodist(exinfo(i), 1);
            exinfo(i).isolation_distance_drug = spkiso.isolation_distance;
            disp(['isolation distance was successfully added to row ' num2str(i)])
        catch
            disp(['row' num2str(i) ' remained nan or emp...'])
        end
%     else
%         disp(['row' num2str(i) ' was skipped (already containing data).'])
%         continue
%     end
end

function [spkiso] = isodist(exinfo, drug)
%% compute spike-isolation quality (S/N, ISI, L, Lratio, islation distance)
% INPUT: fname ... directory of ex-file with c1 or c2 unit with recording data

% +++++++++++++++++++++++++++++++++++++++++++++++++

if drug==0
    fname = exinfo.fname;
else
    fname = exinfo.fname_drug;
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
spkiso.c0.energy = [];
spkiso.c1.energy = [];
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
    
    % energy, store amplitude of waves, store waves    
    for k = 1:length(spkt0)
        energy = sqrt(sum(c0wave(k,:).^2))/length(c0wave(k,:));
        spkiso.c0.energy = [spkiso.c0.energy; energy];
        spkiso.c0.norm_waves = [spkiso.c0.norm_waves; c0wave(k,1:l)/energy];
    end
    for k = 1:length(spkt1)
        energy =  sqrt(sum(c1wave(k,:).^2))/length(c1wave(k,:));
        spkiso.c1.energy = [spkiso.c1.energy; energy];
        spkiso.c1.norm_waves = [spkiso.c1.norm_waves; c1wave(k,1:l)/energy];
    end
end

% PCA for the spike waveform --------------------
len_spk0 = size(spkiso.c0.norm_waves,1);
len_spk1 = size(spkiso.c1.norm_waves,1);

X = [spkiso.c0.norm_waves; spkiso.c1.norm_waves];
sigma = (X'*X)/size(X,1);
[U,S,V] = svd(sigma);
pca = X*U(:,1:3);

% c0
spkiso.c0.norm_pca = pca(1:len_spk0,:);

% c1
spkiso.c1.norm_pca = pca(len_spk0+1:end,:);


% The Mahalanobis distance ------------------------
% c0
X0 = [spkiso.c0.energy, spkiso.c0.norm_pca(:,1), spkiso.c0.time];
% X0 = [spkiso.c0.energy, spkiso.c0.pca(:,1)];
X00 = [spkiso.c0.energy, spkiso.c0.norm_pca(:,1:3), spkiso.c0.time];
% X00 = [spkiso.c0.energy, spkiso.c0.pca(:,1), spkiso.c0.pca(:,2)];

% c1
X1 = [spkiso.c1.energy, spkiso.c1.norm_pca(:,1), spkiso.c1.time];
% X1 = [spkiso.c1.energy, spkiso.c1.pca(:,1)];
X11 = [spkiso.c1.energy, spkiso.c1.norm_pca(:,1:3), spkiso.c1.time];
% X11 = [spkiso.c1.energy, spkiso.c1.pca(:,1), spkiso.c1.pca(:,2)];

% mahalanobis for c0
% % within
% spkiso.c0.mahalanobis(:,1) = mahal(X0, X0);
% spkiso.c0.mahalanobis2(:,1) = mahal(X00, X00);

% outside
spkiso.c0.mahalanobis = mahal(X0, X1);
spkiso.c0.mahalanobis2 = mahal(X00, X11);

% % mahalanobis for c1
% % within
% spkiso.c1.mahalanobis(:,1) = mahal(X1, X1);
% spkiso.c1.mahalanobis2(:,1) = mahal(X11, X11);
% 
% % outside
% spkiso.c1.mahalanobis(:,2) = mahal(X1, X0);
% spkiso.c1.mahalanobis2(:,2) = mahal(X11, X00);


% % Lratio -----------------------------------------------------------
% spkiso.L = [0 0];
% for i = 1:len_spk0
%     spkiso.L(1) = spkiso.L(1) + (1 - chi2cdf(spkiso.c0.mahalanobis(i,2),2));
%     spkiso.L(2) = spkiso.L(2) + (1 - chi2cdf(spkiso.c0.mahalanobis2(i,2),3));
% end
% spkiso.Lratio = spkiso.L/len_spk1;

% Isolation distance ----------------------------------------
if len_spk0 > len_spk1
    sorted_md = sort(spkiso.c0.mahalanobis);
    spkiso.isolation_distance(1) = sorted_md(len_spk1);
    sorted_md = sort(spkiso.c0.mahalanobis2);
    spkiso.isolation_distance(2) = sorted_md(len_spk1);
else
    iso_dist = nan(2,10);
    for r = 1:10
        
        % random sampling
        idx = datasample(1:len_spk1, len_spk0, 'Replace', false);
        
        % PCA for the spike waveform --------------------
        X = [spkiso.c0.norm_waves; spkiso.c1.norm_waves(idx,:)];
        sigma = (X'*X)/size(X,1);
        [U,S,V] = svd(sigma);
        pca = X*U(:,1:2);

        % c0
        norm_pca0 = pca(1:len_spk0,:);

        % c1
        norm_pca1 = pca(len_spk0+1:end,:);

        % The Mahalanobis distance ------------------------
        % c0
        X0 = [spkiso.c0.energy, norm_pca0(:,1)];
        X00 = [spkiso.c0.energy, norm_pca0(:,1), norm_pca0(:,2)];

        % c1
        X1 = [spkiso.c1.energy(idx), norm_pca1(:,1)];
        X11 = [spkiso.c1.energy(idx), norm_pca1(:,1), norm_pca1(:,2)];

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


