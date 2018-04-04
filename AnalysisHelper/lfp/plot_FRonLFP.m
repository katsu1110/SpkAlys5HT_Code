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
lens = length(fr_lfp.session);
for i = 1:lens
    % stLFP in each unit
    figure(1);
    subplot(3,4,i)
    for k = 1:length(fr_lfp.session(i).results.stlfp.vals)-1
        fill_between(
    


