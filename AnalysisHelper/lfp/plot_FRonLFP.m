function plot_FRonLFP
% check whether firing rate modulated by contrast stimulus affects stLFP
% ampltude
%
% Test whether the magnitude of the spike triggered LFP is depending on the
% spiking activity. We use the data recorded with a 2s stimulus, to avoid
% dominant slow fluctuations.
%
% 
% 04.04.18 Katsuhisa wrote it

close all

%% folder specifications
if mean(ismember('gpfs0', cd))==1
    main_dir = '/gpfs01/nienborg/group/Katsuhisa/serotonin_project/LFP_project/MagstLFP_vs_contrast/'; 
else
   main_dir = 'Z:\Katsuhisa\serotonin_project\LFP_project\MagstLFP_vs_contrast\';
end

%% load data given by 'check_FRonLFP.m'
load([main_dir 'Data/fr_lfp.mat'])

%% visualization
% stLFP in each session
lens = length([fr_lfp.session]);
wnd = 0.1; % was 0.3
rho = nan(lens, 2);
c = 0;
for i = 1:lens
    if ~isempty(fr_lfp.session(i).results)
        c = c + 1;
        % stLFP amplitude vs co in each unit
        figure(1);
        subplot(2,5,c)
        plot(fr_lfp.session(i).results.stm.vals(2:end), ...
            fr_lfp.session(i).results.stlfp.peak_stlfp(2:end), ...
            ':ok')
        rho(i,1) = corr(log(fr_lfp.session(i).results.stm.vals(2:end))',...
            fr_lfp.session(i).results.stlfp.peak_stlfp(2:end)', 'type', 'Spearman');
        if ismember(i, 5)
            ylabel('peak stLFP amplitude (uV)')
        end
        if ismember(i, 10)
            xlabel('contrast')
        end
        title(fr_lfp.session(i).results.ntr)
        set(gca, 'XScale', 'log')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

        % stLFP peak time vs co in each unit
        figure(2);
        subplot(2,5,c)
        plot(fr_lfp.session(i).results.stm.vals(2:end), ...
            fr_lfp.session(i).results.stlfp.t_peak_stlfp(2:end), ...
            ':ok')
        rho(i,2) = corr(log(fr_lfp.session(i).results.stm.vals(2:end))',...
            fr_lfp.session(i).results.stlfp.t_peak_stlfp(2:end)', 'type', 'Spearman');
        if ismember(i, 5)
            ylabel('time of peak stLFP (s)')
        end
        if ismember(i, 10)
            xlabel('contrast')
        end
        title(fr_lfp.session(i).results.ntr)
        set(gca, 'XScale', 'log')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

        % fr vs co in each unit
        figure(3);
        subplot(2,5,c)
        plot(fr_lfp.session(i).results.stm.vals(2:end), ...
            fr_lfp.session(i).results.stm.fr(2:end), ...
            ':ok')
        rho(i,2) = corr(log(fr_lfp.session(i).results.stm.vals(2:end))',...
            fr_lfp.session(i).results.stm.fr(2:end)', 'type', 'Spearman');
        if ismember(i, 5)
            ylabel('time of peak stLFP (s)')
        end
        if ismember(i, 10)
            xlabel('contrast')
        end
        title(fr_lfp.session(i).results.ntr)
        set(gca, 'XScale', 'log')
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

        % LFP & stLFP traces
        lenstm = length(fr_lfp.session(i).results.stm.vals(2:end)); 
        col = jet(lenstm);
        pl = cell(1,2);
        leg = cell(1, lenstm);
        for k = 1:lenstm
            figure(4);
            subplot(2,5,c)
            pl{1}(k) = plot(fr_lfp.session(i).results.ts, ...
            fr_lfp.session(i).results.lfp_stm.mean(k+1,:), 'color', col(k,:));
            hold on;
            figure(5);
            subplot(2,5,c)
            fill_between(-wnd:0.001:wnd, ...
                fr_lfp.session(i).results.stlfp.avg_stlfp(k+1, :) - fr_lfp.session(i).results.stlfp.sem_stlfp(k+1, :), ...
                fr_lfp.session(i).results.stlfp.avg_stlfp(k+1, :) + fr_lfp.session(i).results.stlfp.sem_stlfp(k+1, :), ...
                col(k,:), 0.4)
            hold on;
            pl{2}(k) = plot(-wnd:0.001:wnd, fr_lfp.session(i).results.stlfp.avg_stlfp(k+1, :), '-', 'color', col(k,:));
            leg{k} = ['co: ' num2str(fr_lfp.session(i).results.stm.vals(k+1))];
            hold on;
        end
        for ff = 4:5
            figure(ff);
            subplot(2,5,c)
            legend(pl{ff-3}, leg, 'location', 'southeast')
            legend('boxoff')            
            title(fr_lfp.session(i).results.ntr)
            set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
            yy = get(gca, 'YLim');
            plot([0 0], yy, ':k')
            if ff==4
                xlim([fr_lfp.session(i).results.ts(1) fr_lfp.session(i).results.ts(end)])
            elseif ff==5
                xlim([-wnd wnd])
            end
        end
        
        if ismember(i, 5)
            ylabel('stLFP (uV)')
        end
        if ismember(i, 10)
            xlabel('contrast')
        end
    end
end
    
% % across units
% figure(4);
% subplot(1,2,1)
% me_hist(rho(:,1))
% xlabel('Spearman rho')
% title('co vs peak stLFP amp.')
% subplot(1,2,2)
% me_hist(rho(:,2))
% xlabel('Spearman rho')
% title('co vs peak stLFP time.')

