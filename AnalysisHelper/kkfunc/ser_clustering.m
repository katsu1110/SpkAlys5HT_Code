function [datamat, expmat] = ser_clustering(exinfo)
%% try out 2d visualization of multidimensional data in the 5HT experiment

% color-code
cmap = colormap(lines(10));
close all;

% select experiments
load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stim_cond_2007.mat')
load('Z:\Corinna\SharedCode\Katsu\listpvalue_modulation.mat')
exinfo = exinfo(incl_i);
exinfo = exinfo([exinfo.is5HT]==1);

% obtain values from experiments
len_exp = length(exinfo);
datamat = nan(len_exp, 9);
expmat = nan(len_exp, 4);
for i = 1:len_exp
    
    % nonparam_ratio
    datamat(i,1) = exinfo(i).nonparam_ratio;
    
    try
%         if isempty(strfind(exinfo(i).fname, 'all.grating')) && ...
%                                 isempty(strfind(exinfo(i).fname_drug, 'all.grating'))
            % fano factor
            ff_base = exinfo(i).ff_re(1, 1:end-1);
            ff_drug = exinfo(i).ff_re_drug(1, 1:end-1);
            idx = intersect(find(~isnan(ff_base)), find(~isnan(ff_drug)));
            datamat(i,2) = mean(ff_base(idx)) - mean(ff_drug(idx));

            % noise correlation
            datamat(i,3) = (exinfo(i).nc_re(1,1) - 0.0062698*[exinfo(i).c0geomn]) - ...
                (exinfo(i).nc_re_drug(1,1) - 0.0062698*[exinfo(i).c0geomn_drug]);

            % signal correlation
            datamat(i,4) = exinfo(i).rsig - exinfo(i).rsig_drug;
            
            expmat(i,4) = 0;
%         else
%             expmat(i,4) = 1;
%         end
    catch
        expmat(i,4) = 1;
        continue
    end

    % phase selectivity
    datamat(i,5) = mean([exinfo(i).phasesel,exinfo(i).phasesel_drug]);
    
    % RF size
    datamat(i,6) = exinfo(i).RFwx(1);
    datamat(i,7) = exinfo(i).RFwy(1);
    
    % baseline FR
    datamat(i,8) = median(exinfo(i).ratemn);
    
    % waveform width
    datamat(i,9) = mean(exinfo(i).wdt);
    
    % significant unit or not
    expmat(i,1) = p_modulation{i}(2);
    
    % animal
    expmat(i,2) = exinfo(i).ismango;
    
    % stimulus
    switch exinfo(i).param1
        case 'or'
            expmat(i,3) = 1;
        case 'sf'
            expmat(i,3) = 2;
        case 'co'
            expmat(i,3) = 3;
        case 'sz'
            expmat(i,3) = 4;
    end
    
    % contain nans?
    if ismember(1, isnan(datamat(i,:)))
        expmat(i,4) = 1;
    end
end

% visualize via tsne and PCA
figure;
sz = 20;

datamat(expmat(:,4)==1,:) = [];
expmat(expmat(:,4)==1,:) = [];

for v = 1:2
    switch v
        case 1
            subplot(2,2,1)
            % TSNE
            Y = tsne(datamat, 'Algorithm','exact');
            label = 'tsne';
        case 2
            subplot(2,2,2)
            % PCA
            zdatamat = zscore(datamat);
            C = (zdatamat'*zdatamat)/size(zdatamat,2);
            [U,S,V] = svd(C);
            Y = zdatamat*U(:,1:2);
            label = 'pca';
            eigenvals = sum(S,1); 
            varianceExplained = sum(eigenvals(1:2))/sum(eigenvals);
            disp(['variance explained by using 2 PCs: ' num2str(varianceExplained)])
    end

    % loop for visualization based on stimulus types
    for i = 1:4

        % stimulus
        s = expmat(:,3)==i;
        switch i
            case 1
                col = 'r';
            case 2
                col = 'b';
            case 3
                col = 'g';
            case 4
                col = 'm';
        end

        % mango
        Ym = Y(expmat(:,2)==1 & s,:);
        p = expmat(expmat(:,2)==1 & s, 1);

        % significant one
        scatter(Ym(p < 0.05,1), Ym(p < 0.05,2), sz, 'filled', 'marker','o','markerfacecolor',col,'markeredgecolor',col,...
            'markerfacealpha',0.4,'markeredgealpha',0.8)
        hold on;

        % non-significant one
        scatter(Ym(p >= 0.05,1), Ym(p >= 0.05,2), sz, 'filled', 'marker','o','markerfacecolor','w','markeredgecolor',col,...
            'markerfacealpha',0.4,'markeredgealpha',0.8)
        hold on;

        % kaki
        Yk = Y(expmat(:,2)==0 & s,:);
        p = expmat(expmat(:,2)==0 & s, 1);

        % significant one
        scatter(Yk(p < 0.05,1), Yk(p < 0.05,2), sz, 'filled', 'marker','s','markerfacecolor',col,'markeredgecolor',col,...
            'markerfacealpha',0.4,'markeredgealpha',0.8)
        hold on;

        % non-significant one
        scatter(Yk(p >= 0.05,1), Yk(p >= 0.05,2), sz, 'filled', 'marker','s','markerfacecolor','w','markeredgecolor',col,...
            'markerfacealpha',0.4,'markeredgealpha',0.8)
    end

    title(label)
    set(gca,'box','off')
    set(gca,'TickDir','out')
%     axis square
    
    % k-means clustering
    kn = 2;
    idx = [];
    scand = 0;
    disp([label ' ---------------'])
    for k = 2:7
        cidx = kmeans(Y,k,'distance','cityblock','MaxIter',1000);
        silh = mean(silhouette(Y,cidx,'cityblock'));
        disp(['k = ' num2str(k) ', silh = ' num2str(silh)])
        if scand < silh
            scand = silh;
            kn = k;
            idx = cidx;
        end
    end
    
    subplot(2,2,v+2)
    for kk = 1:kn        
        plot(Y(idx==kk, 1), Y(idx==kk, 2), '.', 'color', cmap(kk,:), 'markersize',12);
        hold on;
    end
    title(['the number of cluster = ' num2str(kn) ', silh = ' num2str(scand)])    
    set(gca,'box','off')
    set(gca,'TickDir','out')
%     axis square
end
% set(gcf,'position',[598   411   341   226])