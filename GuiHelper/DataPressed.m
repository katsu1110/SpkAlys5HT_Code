function DataPressed(~, ~, datinfo, fig2plot)
% DataPressed(source, event, datinfo, fig2plot)
% 
% 
% Callback function that opens the figures specified as string in the cell
% fig2plot for the pressed data. 
%
% This function can also be called as a normal function with empty source
% and event input.
% Example:
% DataPressed([], [], exinfo(270), {'LFP Gui'})
% 
%
% Figure names:
% -  'LFP Gui'  opens the LFP Gui figures and spectrograms 
% -  'Variability' opens the figure with z-scored trial spike count and/or
%       noise correlation plot (only DG)
% - 'Tuning Curve' opens the figure with tuning curves
% - 'Wave Form' shows the average spike waveform in each condition
% - 'Regression' opens the type-II regression fit between baseline and drug
% - 'ISI' shows the ACF
% - 'Raster' opens the raster plot (only DG)
% - 'Phase Select.' shows the psth superimposed with the temporal frequency
%   of the stimulus and the computed f1/f0. It shows the results of the
%   experiment with TF if there was one. Otherwise it shows the results of
%   this particular experiment.
% - 'smooth PSTH' shows the smoothed PSTH (only DG)
% - 'Spike Density' opens the figure with the spike density function and
% the deviation of orientation selective components used to estimate the
% response latency (only Reverse Correlation Analysis data)
% -  'Recovery' ?
%
% @CL
% last modified 17.05.2017



disp(datinfo(1).fname_drug);
fprintf('Anova Baseline p=%1.3f, Drug: p=%1.3f \n', datinfo(1).p_anova, datinfo(1).p_anova_drug);
fprintf('recovery overlap p=%1.3f \n', datinfo(1).ret2base);
if strcmp(datinfo(1).param1, 'co')
    fprintf('Add Anova Baseline data<=c50 p=%1.3f, data>c50: p=%1.3f \n', datinfo(1).fitparam.undersmpl(1), datinfo(1).fitparam.undersmpl(2));
    fprintf('Add Anova Drug data<=c50 p=%1.3f, data>c50: p=%1.3f \n', datinfo(1).fitparam_drug.undersmpl(1), datinfo(1).fitparam.undersmpl(2));
end

disp(datinfo(1).spkqual_base);
disp(datinfo(1).spkqual_drug);



j = 1;
pos = 1;
nfig = length(fig2plot);

if any(strcmp(fig2plot, 'LFP Gui'))
    nfig = nfig-1;
end
if any(strcmp(fig2plot, 'Raster'))
    nfig = nfig-1;
end
if any(strcmp(fig2plot, 'Variability'))
    nfig = nfig-1;
end

if any(strcmp(fig2plot, 'Regression'))
    nfig = nfig-1;
end


if nfig > 0
    h = figure('Visible', 'off');
end


temp =[];
while j<=length(fig2plot)
    
    switch fig2plot{j}
       
        case 'Direction Selectivity'
           
           temp = figure;
           getDIRsel( datinfo(1).fname, 'plot')
           set(findobj(gca, 'type', 'line'), 'LineStyle', '-');
           getDIRsel( datinfo(1).fname_drug, 'plot')
           title(sprintf(['DS baseline: %1.2f \n', datinfo(1).drugname ': %1.2f'], ...
               datinfo(1).ds2, datinfo(1).ds2));
                  
           set(findobj(gca, 'type', 'line'), 'Color', getCol(datinfo(1)));
           
        case 'LFP Gui'
            
            ex1 = loadCluster(datinfo(1).fname, 'ocul', datinfo(1).ocul);
            ex2 = loadCluster(datinfo(1).fname_drug, 'ocul', datinfo(1).ocul);
           
            hlfp = LFPGui( ex1, ex2, datinfo(1).figname);
            j = j+1;
            continue
            
        case 'Variability'
            
            col = lines(3);
            figure('UserData', datinfo(1));
            ff_name = 'ff';
            
            subplot(1,2,1);
            scatter([datinfo(1).(ff_name).classic.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic.spkcnt_var], 40, 'filled', ...
                'MarkerFaceColor', col(1,:),...
                'MarkerFaceAlpha', 0.6); hold on;
            
            text([datinfo(1).(ff_name).classic.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic.spkcnt_var], ...
                num2str([datinfo(1).(ff_name).classic.stimrep]'),...
                'Color', col(1,:));
            
            
            scatter([datinfo(1).(ff_name).classic_2ndhalf.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic_2ndhalf.spkcnt_var], 40, 'filled', ...
                'MarkerFaceColor', col(2,:),...
                'MarkerFaceAlpha', 0.6); hold on;
            
            text([datinfo(1).(ff_name).classic_2ndhalf.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic_2ndhalf.spkcnt_var], ...
                num2str([datinfo(1).(ff_name).classic_2ndhalf.stimrep]'),...
                'Color', col(2,:));

            
            scatter([datinfo(1).(ff_name).classic_20plus.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic_20plus.spkcnt_var], 40, 'filled', ...
                'MarkerFaceColor', col(3,:),...
                'MarkerFaceAlpha', 0.6); hold on;
            
            text([datinfo(1).(ff_name).classic_20plus.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic_20plus.spkcnt_var], ...
                num2str([datinfo(1).(ff_name).classic_20plus.stimrep]'), ...
                'Color', col(3,:));
            
            eqax; unity
            xlabel('spike count mean') ;ylabel('spike count variance');
            axis square
            xlim_= get(gca, 'XLim');
            set(gca, 'XLim', [0.1 xlim_(2)]);
            set(gca, 'YLim', [0.1 xlim_(2)]);
            
            subplot(1,2,2);
            ff_name = 'ff_drug';
            scatter([datinfo(1).(ff_name).classic.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic.spkcnt_var], 40, 'filled', ...
                'MarkerFaceColor', col(1,:),...
                'MarkerFaceAlpha', 0.6); hold on;
            
            text([datinfo(1).(ff_name).classic.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic.spkcnt_var], ...
                num2str([datinfo(1).(ff_name).classic.stimrep]'),...
                'Color', col(1,:));
            
            
            scatter([datinfo(1).(ff_name).classic_2ndhalf.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic_2ndhalf.spkcnt_var], 40, 'filled', ...
                'MarkerFaceColor', col(2,:),...
                'MarkerFaceAlpha', 0.6); hold on;
            
            text([datinfo(1).(ff_name).classic_2ndhalf.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic_2ndhalf.spkcnt_var], ...
                num2str([datinfo(1).(ff_name).classic_2ndhalf.stimrep]'),...
                'Color', col(2,:));

            
            scatter([datinfo(1).(ff_name).classic_20plus.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic_20plus.spkcnt_var], 40, 'filled', ...
                'MarkerFaceColor', col(3,:),...
                'MarkerFaceAlpha', 0.6); hold on;
            
            text([datinfo(1).(ff_name).classic_20plus.spkcnt_mn], ...
                [datinfo(1).(ff_name).classic_20plus.spkcnt_var], ...
                num2str([datinfo(1).(ff_name).classic_20plus.stimrep]'), ...
                'Color', col(3,:));
            
            xlabel('spike count mean') ;ylabel('spike count variance');
            axis square
            eqax; unity
            xlim_= get(gca, 'XLim');
            set(gca, 'XLim', [0.1 xlim_(2)]);
            set(gca, 'YLim', [0.1 xlim_(2)]);

            legend('', 'full exp', '2nd half', '20+');

            
            h_var = openfig(datinfo(1).fig_varxtime, 'visible');
            set(h_var, 'UserData', datinfo(1));
            
%             h_var2 = openfig(strrep(datinfo(1).fig_varxtime, 'Variability', 'Phase'), 'visible');
            j = j+1;
            pos = pos+1;
            continue
            
           
        case 'Tuning Curve'
        temp = openfig(datinfo(1).fig_tc, 'invisible');
            
        case 'Wave Form'
            temp = openfig(datinfo(1).fig_waveform, 'invisible');
            
        case 'Regression'
%             temp = openfig(datinfo(1).fig_regl, 'invisible');
            openfig(datinfo(1).fig_regl);
             j = j+1;
            pos = pos+1;
            continue
            
        case 'pval mod'
            temp = openfig(datinfo(1).fig_bri, 'invisible');
            
        case 'Raster'
            
            h_raster = openfig(datinfo(1).fig_raster);
            set(h_raster, 'UserData', datinfo(1));
            j = j+1;

            continue
            
        case 'Phase Select.'
            temp = openfig(datinfo(1).fig_phase);  
            if exist( datinfo(1).fig_phasetf )
                temp = openfig(datinfo(1).fig_phasetf);
            end
            j = j+1; 

            continue
        case 'smooth PSTH'
            temp = openfig(datinfo(1).fig_psth);
            ax = findobj(temp, 'Type', 'Axes');
            delete(ax([2,4]))
            
            ylim_ = [0 max( horzcat(ax([3,5]).YLim) )];
            set(ax([3,5]), 'YLim', ylim_);
                        
        case 'Spike Density'
            temp = openfig(datinfo(1).fig_sdfs, 'invisible');
            
        case 'Recovery'
            
            temp = openfig(datinfo(1).fig_recovery);
            j = j+1;
            continue
           
    end
    
    
    
    ax = findobj(temp, 'Type', 'axes');
    newax = copyobj(ax, h);    
    
    close(temp);
    
    for i = 1:length(newax)
        newax(i).Position = [0.1+(0.9/nfig)*(pos-1)  ...
            0.12+(0.8/length(newax))*(i-1) ...
            0.85/nfig-0.05 ...
            0.8/length(newax)-0.1];
        newax(i).Title.FontSize = 8;
    end
    
    pos = pos+1;
    j = j+1;
end


if nfig>0
    h.Position = [287   500   350*nfig  450];
    h.Name = datinfo(1).figname;
    h.UserData = datinfo;
    h.Visible = 'on';
end
end