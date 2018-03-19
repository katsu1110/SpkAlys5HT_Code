function plotLFP(LFPinfo, exinfo, stmtype)
%% plot LFP data across sessions
%
% load('Z:\Katsuhisa\LFP_project\Data\LFPinfo.mat')
% load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\dataset\Data\exinfo.mat')
%
% written by Katsuhisa (19.03.18)
% ++++++++++++++++++++++++++++++++

close all;

%% LFP data: baseline vs drug
% stimulus type
row = [];
is5ht = [];
if strcmp(stmtype, 'rc')
    load('Z:\Corinna\SharedCode\Katsu\list_RC.mat')
    for i = 1:length(list_RC)
        if LFPinfo.session(list_RC(i)).exist==1
            row = [row, list_RC(i)];
            if strcmp(exinfo(list_RC(i)).drugname, '5HT')
                is5ht = [is5ht, 1];
            else
                is5ht = [is5ht, 0];
            end
        end
    end
elseif strcmp(stmtype, 'co')
    load('Z:\Corinna\SharedCode\Katsu\incl_i_all_stim_cond_2007.mat')
    for i = 1:length(incl_i)
        if strcmp(exinfo(incl_i(i)).param1, stmtype)
            row = [row, incl_i(i)];
            if strcmp(exinfo(incl_i(i)).drugname, '5HT')
                is5ht = [is5ht, 1];
            else
                is5ht = [is5ht, 0];
            end
        end
    end
end

lenr = length(row);
    
for f = 1:3
    switch f
        case 1
            fieldname = 'drug';
        case 2
            fieldname = 'ps_base';
        case 3
            fieldname = 'ps_drug';
    end
    
    % stLFP ==========================
    % LFP
    h = figure;
    len_trace = length(LFPinfo.session(row(1)).results.(fieldname).cond(1).stlfp.mean);
    len_pow = length(LFPinfo.session(row(1)).results.(fieldname).cond(1).stlfp.pow);
    for k = 1:2
        para.cond(k).stlfp.zeroval = nan(1, lenr);
        para.cond(k).stlfp.trace = nan(lenr, len_trace);
        para.cond(k).stlfp.power = nan(lenr, len_pow);
    end
    for i = 1:lenr
        for k = 1:2
            para.cond(k).stlfp.zeroval(i) = ...
                LFPinfo.session(row(i)).results.(fieldname).cond(k).stlfp.mean(round(len_trace/2));
            para.cond(k).stlfp.trace(i,:) = ...
                LFPinfo.session(row(i)).results.(fieldname).cond(k).stlfp.mean;
            para.cond(k).stlfp.pow(i,:) = ...
                LFPinfo.session(row(i)).results.(fieldname).cond(k).stlfp.pow';
        end
    end
    for k = 1:2
        switch k
            case 1
                drugname = 'NaCl';
            case 2
                drugname = '5HT';
        end

        % stLFP traces --- 5HT
        subplot(2,3,1+3*(k-1))
        wnd = 0.064;
        me = nanmean(para.cond(1).stlfp.trace(is5ht==k-1,:), 1);
        sem = nanstd(para.cond(1).stlfp.trace(is5ht==k-1,:), [], 1)...
            /sqrt(sum(is5ht==k-1));
        fill_between(-wnd:0.001:wnd, me - sem, me + sem, zeros(1,3))
        hold on;
        plot(-wnd:0.001:wnd, me, '-k')
        hold on;
        me = nanmean(para.cond(2).stlfp.trace(is5ht==k-1,:), 1);
        sem = nanstd(para.cond(2).stlfp.trace(is5ht==k-1,:), [], 1)...
            /sqrt(sum(is5ht==k-1));
        fill_between(-wnd:0.001:wnd, me - sem, me + sem, [1 0 0])
        hold on;
        plot(-wnd:0.001:wnd, me, '-r')
        hold on;
        yy = get(gca, 'YLim');
        plot([0 0],yy, '-k')
        xlim([-wnd wnd])
        ylim(yy)
        xlabel('time (s)')
        ylabel('LFP')
        title(drugname)
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square

        % stLFP power
        subplot(2,3,2+3*(k-1))
        freq = LFPinfo.session(row(1)).results.(fieldname).cond(1).stlfp.freq;
        me = nanmean(para.cond(1).stlfp.pow(is5ht==k-1,:), 1);
        sem = nanstd(para.cond(1).stlfp.pow(is5ht==k-1,:), [], 1)...
            /sqrt(sum(is5ht==k-1));
        fill_between(freq, me - sem, me + sem, zeros(1,3))
        hold on;
        plot(freq, me, '-k')
        hold on;
        me = nanmean(para.cond(2).stlfp.pow(is5ht==k-1,:), 1);
        sem = nanstd(para.cond(2).stlfp.pow(is5ht==k-1,:), [], 1)...
            /sqrt(sum(is5ht==k-1));
        fill_between(freq, me - sem, me + sem, [1 0 0])
        hold on;
        plot(freq, me, '-r')
        hold on;
        xlim([min(freq) max(freq)])
        xlabel('frequency (Hz)')
        ylabel('power')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square

        % stLFP scatter
        subplot(2,3,3)
        unity_plot(para.cond(1).stlfp.zeroval, para.cond(2).stlfp.zeroval, is5ht)
        title('stLFP')
        xlabel('baseline')
        ylabel('drug')
    end
    set(h, 'Name', [fieldname ': stLFP'], 'NumberTitle','off')

    % LFP power
    h = figure;
    for k = 1:2
        para.cond(k).stlfp.power_delta = nan(1, lenr);
        para.cond(k).stlfp.power_theta = nan(1, lenr);
        para.cond(k).stlfp.power_alpha = nan(1, lenr);
        para.cond(k).stlfp.power_beta = nan(1, lenr);
        para.cond(k).stlfp.power_gamma = nan(1, lenr);
    end
    for i = 1:lenr
        for k = 1:2
            pow = LFPinfo.session(row(i)).results.(fieldname).cond(k).stlfp.pow;
            freq = LFPinfo.session(row(i)).results.(fieldname).cond(k).stlfp.freq;
            para.cond(k).stlfp.power_delta(i) = mean(pow(freq > 0 & freq < 4));
            para.cond(k).stlfp.power_theta(i) = mean(pow(freq >= 4 & freq < 8));
            para.cond(k).stlfp.power_alpha(i) = mean(pow(freq >= 8 & freq <= 13));
            para.cond(k).stlfp.power_beta(i) = mean(pow(freq >= 14 & freq <= 29));
            para.cond(k).stlfp.power_gamma(i) = mean(pow(freq >= 30 & freq <= 80));
        end
    end
    for b = 1:5
        switch b
            case 1
                band = 'delta';
            case 2
                band = 'theta';
            case 3
                band = 'alpha';
            case 4
                band = 'beta';
            case 5
                band = 'gamma';
        end
        subplot(2,3,b)
        unity_plot(para.cond(1).stlfp.(['power_' band]),...
            para.cond(2).stlfp.(['power_' band]), is5ht)
        title(band)
    end
    set(h, 'Name', [fieldname ': stLFP power bands'], 'NumberTitle','off')

    % LFP tuning to stimuli =================
    h = figure;
    for k = 1:2
        para.cond(k).lfpstm.power_delta = nan(1, lenr);
        para.cond(k).lfpstm.power_theta = nan(1, lenr);
        para.cond(k).lfpstm.power_alpha = nan(1, lenr);
        para.cond(k).lfpstm.power_beta = nan(1, lenr);
        para.cond(k).lfpstm.power_gamma = nan(1, lenr);
    end
    stm = [];
    for i = 1:lenr
        stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
        stm = [stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx)];
        for k = 1:2       
            for b = 1:5
                switch b
                    case 1
                        band = 'delta';
                    case 2
                        band = 'theta';
                    case 3
                        band = 'alpha';
                    case 4
                        band = 'beta';
                    case 5
                        band = 'gamma';
                end
                para.cond(k).lfpstm.power.(band)(i) = ...
                    mean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.lfp_stm_wave(b).power(stmidx));
            end
        end
    end
    [stm, counts] = uniquecount(stm);
    stm(counts < mean(counts)*0.3) = [];
    lenstm = length(stm);
    for k = 1:2
        para.cond(k).lfpstm.power_tuning.delta = nan(lenr, lenstm);
        para.cond(k).lfpstm.power_tuning.theta = nan(lenr, lenstm);
        para.cond(k).lfpstm.power_tuning.alpha = nan(lenr, lenstm);
        para.cond(k).lfpstm.power_tuning.beta = nan(lenr, lenstm);
        para.cond(k).lfpstm.power_tuning.gamma = nan(lenr, lenstm);
    end
    for i = 1:lenr
        stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
        idx = ismember(LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx), stm);
        for b = 1:5
            switch b
                case 1
                    band = 'delta';
                case 2
                    band = 'theta';
                case 3
                    band = 'alpha';
                case 4
                    band = 'beta';
                case 5
                    band = 'gamma';
            end

            pow = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.lfp_stm_wave(b).power(stmidx);
            para.cond(1).lfpstm.power_tuning.(band)(i,:) = pow(idx);
            pow = LFPinfo.session(row(i)).results.(fieldname).cond(2).lfpstm.lfp_stm_wave(b).power(stmidx);
            para.cond(2).lfpstm.power_tuning.(band)(i,:) = pow(idx);

            % normalization
            m = mean([para.cond(1).lfpstm.power_tuning.(band)(i,:), ...
                para.cond(2).lfpstm.power_tuning.(band)(i,:)]);
            para.cond(1).lfpstm.power_tuning.(band)(i,:) = ...
                para.cond(1).lfpstm.power_tuning.(band)(i,:)/m;
            para.cond(2).lfpstm.power_tuning.(band)(i,:) = ...
                para.cond(2).lfpstm.power_tuning.(band)(i,:)/m;
        end
    end
    for b = 1:5
        switch b
            case 1
                band = 'delta';
            case 2
                band = 'theta';
            case 3
                band = 'alpha';
            case 4
                band = 'beta';
            case 5
                band = 'gamma';
        end
        % effect size of 5HT on tuning
        subplot(2,5,b)
        unity_plot(para.cond(1).lfpstm.power.(band),...
            para.cond(2).lfpstm.power.(band), is5ht)
        title(band)
        if b==1
            xlabel('baseline')
            ylabel('drug')
        end

        % tuning curve
        for l = 1:2
            switch l
                case 1
                    col = [0 0 0];
                case 2
                    col = [1 0 0];
            end
            subplot(2,5,b+5)
            me = mean(para.cond(1).lfpstm.power_tuning.(band)(is5ht==l-1,:), 1);
            sem = std(para.cond(1).lfpstm.power_tuning.(band)(is5ht==l-1,:), [], 1)/...
                sqrt(size(para.cond(1).lfpstm.power_tuning.(band)(is5ht==l-1,:),1));
            errorbar(stm, me, sem, 'color', col, 'linestyle','-', 'CapSize', 0)
            hold on;
            me = mean(para.cond(2).lfpstm.power_tuning.(band)(is5ht==l-1,:), 1);
            sem = std(para.cond(2).lfpstm.power_tuning.(band)(is5ht==l-1,:), [], 1)/...
                sqrt(size(para.cond(2).lfpstm.power_tuning.(band)(is5ht==l-1,:),1));
            errorbar(stm, me, sem, 'color', col, 'linestyle','--', 'CapSize', 0)
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
        end
        if b==1
            xlabel(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stm.param)
            ylabel('normalized power')
        end
    end
    set(h, 'Name', [fieldname ': Tuning of LFP power'], 'NumberTitle','off')


    % stLFP tuning to stimuli =================
    h = figure;
    for k = 1:2
        para.cond(k).lfpstm.stlfppower.delta = nan(1, lenr);
        para.cond(k).lfpstm.stlfppower.theta = nan(1, lenr);
        para.cond(k).lfpstm.stlfppower.alpha = nan(1, lenr);
        para.cond(k).lfpstm.stlfppower.beta = nan(1, lenr);
        para.cond(k).lfpstm.stlfppower.gamma = nan(1, lenr);
    end
    stm = [];
    for i = 1:lenr
        stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
        stm = [stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx)];

        for k = 1:2    
            para.cond(k).lfpstm.stlfppower.delta(i) = ...
                nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,1), 1);
            para.cond(k).lfpstm.stlfppower.theta(i) = ...
                nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,2), 1);
            para.cond(k).lfpstm.stlfppower.alpha(i) = ...
                nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,3), 1);
            para.cond(k).lfpstm.stlfppower.beta(i) = ...
                nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,4), 1);
            para.cond(k).lfpstm.stlfppower.gamma(i) = ...
                nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,5), 1);
        end
    end
    [stm, counts] = uniquecount(stm);
    stm(counts < mean(counts)*0.3) = [];
    lenstm = length(stm);
    for k = 1:2
        para.cond(k).lfpstm.stlfppower_tuning.delta = nan(lenr, lenstm);
        para.cond(k).lfpstm.stlfppower_tuning.theta = nan(lenr, lenstm);
        para.cond(k).lfpstm.stlfppower_tuning.alpha = nan(lenr, lenstm);
        para.cond(k).lfpstm.stlfppower_tuning.beta = nan(lenr, lenstm);
        para.cond(k).lfpstm.stlfppower_tuning.gamma = nan(lenr, lenstm);
    end
    for i = 1:lenr
        stmidx = find(LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000);
        idx = ismember(LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx), stm);
        idx2 = ismember(stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx));

        for k = 1:2    
            bandtu = LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.stlfp.band(stmidx,:);
            para.cond(k).lfpstm.stlfppower_tuning.delta(i,idx2) = bandtu(idx, 1)';
            para.cond(k).lfpstm.stlfppower_tuning.theta(i,idx2) = bandtu(idx, 2)';
            para.cond(k).lfpstm.stlfppower_tuning.alpha(i,idx2) = bandtu(idx, 3)';
            para.cond(k).lfpstm.stlfppower_tuning.beta(i,idx2) = bandtu(idx, 4)';
            para.cond(k).lfpstm.stlfppower_tuning.gamma(i,idx2) = bandtu(idx, 5)';
        end

        for b = 1:5
            switch b
                case 1
                    band = 'delta';
                case 2
                    band = 'theta';
                case 3
                    band = 'alpha';
                case 4
                    band = 'beta';
                case 5
                    band = 'gamma';
           end

            % normalization
            m = nanmean([para.cond(1).lfpstm.stlfppower_tuning.(band)(i,:), ...
                para.cond(2).lfpstm.stlfppower_tuning.(band)(i,:)]);
            para.cond(1).lfpstm.stlfppower_tuning.(band)(i,:) = ...
                para.cond(1).lfpstm.stlfppower_tuning.(band)(i,:)/m;
            para.cond(2).lfpstm.stlfppower_tuning.(band)(i,:) = ...
                para.cond(2).lfpstm.stlfppower_tuning.(band)(i,:)/m;
        end
    end

    for b = 1:5
        switch b
            case 1
                band = 'delta';
            case 2
                band = 'theta';
            case 3
                band = 'alpha';
            case 4
                band = 'beta';
            case 5
                band = 'gamma';
        end

        % effect size of 5HT on tuning
        subplot(2,5,b)
        unity_plot(para.cond(1).lfpstm.stlfppower.(band),...
            para.cond(2).lfpstm.stlfppower.(band), is5ht)
        title(band)
        if b==1
            xlabel('baseline')
            ylabel('drug')
        end

        % tuning curve
        for l = 1:2
            switch l
                case 1
                    col = [0 0 0];
                case 2
                    col = [1 0 0];
            end
            subplot(2,5,b+5)
            me =nanmean(para.cond(1).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:), 1);
            sem = nanstd(para.cond(1).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:), [], 1)/...
                sqrt(size(para.cond(1).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:),1));
            errorbar(stm, me, sem, 'color', col, 'linestyle','-', 'CapSize', 0)
            hold on;
            me = nanmean(para.cond(2).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:), 1);
            sem = nanstd(para.cond(2).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:), [], 1)/...
                sqrt(size(para.cond(2).lfpstm.stlfppower_tuning.(band)(is5ht==l-1,:),1));
            errorbar(stm, me, sem, 'color', col, 'linestyle','--', 'CapSize', 0)
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
        end
        if b==1
            xlabel(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stm.param)
            ylabel('normalized power')
        end
    end
    set(h, 'Name', [fieldname ': Tuning of stLFP power'], 'NumberTitle','off')


    % spike-LFP coherence ==================
    h = figure;
    stm = [];
    for i = 1:lenr
        stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
        stm = [stm, LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx)];
    end
    [stm, counts] = uniquecount(stm);
    stm(counts < mean(counts)*0.3) = [];
    lenstm = length(stm);
    freq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.coherence.f{1};
    lenf = length(freq);
    for k = 1:2
        para.cond(k).lfpstm.coherence_stmval = nan(lenr, lenstm);
        para.cond(k).lfpstm.coherence_freq = nan(lenr, lenf);
    end
    for i = 1:lenr
        stmidx = LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals < 1000;
        idx = ismember(LFPinfo.session(row(i)).results.(fieldname).cond(1).lfpstm.stm.vals(stmidx), stm);
        idxi = find(idx==1);
        for k = 1:2
            for s = 1:sum(idx==1)
                para.cond(k).lfpstm.coherence_stmval(i, s) = ...
                    nanmean(LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.coherence.C{idxi(s)});      
            end
            para.cond(k).lfpstm.coherence_freq(i,:) = ...
                    LFPinfo.session(row(i)).results.(fieldname).cond(k).lfpstm.coherence.C{idxi(s)}';
        end
    end

    for l = 1:2
        switch l
            case 1
                drugname = 'NaCl';
            case 2
                drugname = '5HT';
        end

        % coherence vs frequency
        subplot(2,2,1+2*(l-1))
        freq = LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.coherence.f{1};
        me = nanmean(para.cond(1).lfpstm.coherence_freq(is5ht==l-1,:),1);
        sem = nanstd(para.cond(1).lfpstm.coherence_freq(is5ht==l-1,:), [], 1)...
            /sqrt(lenr);
        fill_between(freq, me - sem, me + sem, [0 0 0])
        hold on;
        plot(freq, me, '-k')
        hold on;
        me = nanmean(para.cond(2).lfpstm.coherence_freq(is5ht==l-1,:),1);
        sem = nanstd(para.cond(2).lfpstm.coherence_freq(is5ht==l-1,:), [], 1)...
            /sqrt(lenr);
        fill_between(freq, me - sem, me + sem, [1 0 0])
        hold on;
        plot(freq, me, '-r')
        hold on;
        xlabel('frequency (Hz)')
        ylabel('coherence')
        title(drugname)
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square

        % tuning curve
        switch l
            case 1
                col = [0 0 0];
            case 2
                col = [1 0 0];
        end

        subplot(2,2,2+2*(l-1))
        me = nanmean(para.cond(1).lfpstm.coherence_stmval(is5ht==l-1,:),1);
        sem = nanstd(para.cond(1).lfpstm.coherence_stmval(is5ht==l-1,:), [], 1)/sqrt(lenr);        
        errorbar(stm, me, sem, 'color', col, 'linestyle','-', 'CapSize', 0)
        hold on;
        me = nanmean(para.cond(2).lfpstm.coherence_stmval(is5ht==l-1,:),1);
        sem = nanstd(para.cond(2).lfpstm.coherence_stmval(is5ht==l-1,:), [], 1)/sqrt(lenr);
        errorbar(stm, me, sem, 'color', col, 'linestyle','--', 'CapSize', 0)
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square
        xlabel(LFPinfo.session(row(1)).results.(fieldname).cond(1).lfpstm.stm.param)
    end

    set(h, 'Name', [fieldname ': spike-LFP coherence'], 'NumberTitle','off')
end


%% subfunctions
function unity_plot(x,y,bi)
if nargin < 3
    bi = ones(length(x), 1);
end
oks = ~isnan(x) & ~isnan(y);
x = x(oks);
y = y(oks);
bi = bi(oks);
minima = min([x, y]); 
maxima = max([x, y]); 
dist = maxima - minima;
minima = minima - dist*0.1;
maxima = maxima + dist*0.1;
plot([minima maxima], [minima maxima], '-', 'color', 0.6*ones(1,3))
hold on;
for i = 1:length(bi)
    if bi(i)==1
        col = [1 0 0];
    else
        col = [0 0 0];
    end
    scatter(x(i), y(i), 20, 'o',...
            'markerfacecolor', col, 'markeredgecolor', col, ...
            'markerfacealpha', 0.4, 'markeredgealpha', 0.8);
        hold on;
end
if nargin < 3
    try
        p2 = signrank(x, y);
    catch
        p2 = nan;
    end
    text(minima + 0.7*dist, minima + 0.15*dist, ['p = ' num2str(p2)])
else
    try
        p1 = signrank(x(bi==0),y(bi==0));
        p2 = signrank(x(bi==1),y(bi==1));
    catch
        p1 = nan;
        p2 = nan;
    end
    text(minima + 0.7*dist, minima + 0.15*dist, ['nacl: p = ' num2str(p1)])
    text(minima + 0.7*dist, minima + 0.25*dist, ['5ht: p = ' num2str(p2)])
end

axis([minima maxima minima maxima])
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out'); axis square

