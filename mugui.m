function mugui(exinfo)



%% add the results of the RC experiment subsampling for noise corrected latency estimates
exinfo = addLM2Struct(exinfo); 


%% predefine inclusion criteria
nrep = 4;
minspk = 10;

%% predefine variables

suFigsOnDataClick = {}; % cell array of strings that specify which plots are shown when a single data point is clicked
dat = [];   % structure that contains the data to be plotted and the original exinfo strucure
incl_i = 1:length(exinfo); % indices corresponding to the entries in exinfo that are supposed to be shown in the plot, see UpdateInclusion for more information
margins = 0.03;

%% loading GUI information 
guiprop = PlotProps();

axisopts = guiprop.axisOpts; % axis options, values to plot on the x or y axis
datacrit = guiprop.datacrit; % data criteria such as highest dose. specified in PlotProps.m 
electrodecrit = guiprop.electrodecrit;  % electrode critera
stimulicond = guiprop.stimulicond;  % stimulus conditions 
eyeopts = guiprop.eyeopts;  % ocular condition


%% creating the GUI 

% specify the figure size, position, etc.
sz = get(0, 'Screensize');
fig_h = figure('Position', [sz(3)*0.2 sz(4)*0.25 sz(3)*0.6 sz(4)*0.6], ...
    'CloseRequestFcn', @my_closereq);


pos = get(fig_h, 'position');

% all the buttons, checkboxes etc. are created in the following
%------------------------- popmenu for Y
popY_h = uicontrol(fig_h, ....
    'Style',  'popupmenu',...
    'String', [{'y'} axisopts],...
    'Position', [pos(3)*0.2 pos(4)*0.8 pos(3)*0.15 pos(4)*0.1],...
    'Callback', '');

%-------------------------- popmenu for X
popX_h = uicontrol(fig_h, ....
    'Style',  'popupmenu',...
    'String', [{'x'} axisopts],...
    'Position', [pos(3)*0.2 pos(4)*0.6 pos(3)*0.15 pos(4)*0.1],...
    'Callback', '');

%------------------ checkbox for inclusion
addUnity = uicontrol(fig_h, ....
    'Style',  'checkbox',...
    'String', 'unity',...
    'Position', [pos(3)*margins pos(4)*0.61 pos(3)*0.15 pos(4)*0.1]);

addCross = uicontrol(fig_h, ....
    'Style',  'checkbox',...
    'String', 'cross',...
    'Position', [pos(3)*margins pos(4)*0.55 pos(3)*0.15 pos(4)*0.1]);

addRegress = uicontrol(fig_h, ....
    'Style',  'checkbox',...
    'String', 'regression',...
    'Position', [pos(3)*margins pos(4)*0.49 pos(3)*0.15 pos(4)*0.1]);

eqAxes = uicontrol(fig_h, ....
    'Style',  'checkbox',...
    'String', 'equal axes',...
    'Position', [pos(3)*margins pos(4)*0.43 pos(3)*0.15 pos(4)*0.1]);

addHistograms = uicontrol(fig_h, ....
    'Style',  'checkbox',...
    'String', 'histogram',...
    'Position', [pos(3)*margins pos(4)*0.37 pos(3)*0.15 pos(4)*0.1]);


%------------------------- restriction to isolation quality
txt0 = uicontrol(fig_h, 'Style','text',...
    'Position',[pos(3)*margins pos(4)*0.26 pos(3)*0.15 pos(4)*0.1],...
    'String','SpkSort Q', ...
    'HorizontalAlignment', 'left');
editSpkSortQualtiy_h = uicontrol(fig_h, 'Style','edit',...
    'Position',[pos(3)*(margins+0.06) pos(4)*0.335 pos(3)*0.04 pos(4)*0.03], ...
    'String', ' >= 1 ');

%------------------------- r2 gain fit restriction
txt = uicontrol(fig_h, 'Style','text',...
    'Position',[pos(3)*margins pos(4)*0.21 pos(3)*0.15 pos(4)*0.1],...
    'String','R2 reg', ...
    'HorizontalAlignment', 'left');
editr2_h = uicontrol(fig_h, 'Style','edit',...
    'Position',[pos(3)*(margins+0.06) pos(4)*0.285 pos(3)*0.04 pos(4)*0.03], ...
    'String', ' > -inf ');

%------------------------- gain restriction
txt2 = uicontrol(fig_h, 'Style','text',...
    'Position',[pos(3)*margins pos(4)*0.16 pos(3)*0.15 pos(4)*0.1],...
    'String','latency ', ...
    'HorizontalAlignment', 'left');
editlatency_h = uicontrol(fig_h, 'Style','edit',...
    'Position',[pos(3)*(margins+0.06) pos(4)*0.235 pos(3)*0.04 pos(4)*0.03], ...
    'String', ' > -inf ');


%------------------------- r2 gauss fit restriction
txt3 = uicontrol(fig_h, 'Style','text',...
    'Position',[pos(3)*margins pos(4)*0.11 pos(3)*0.15 pos(4)*0.1],...
    'String','R2 TC ', ...
    'HorizontalAlignment', 'left');
editr2tc_h = uicontrol(fig_h, 'Style','edit',...
    'Position',[pos(3)*(margins+0.06) pos(4)*0.185 pos(3)*0.04 pos(4)*0.03], ...
    'String', ' > -inf ');

%------------------------- tc anova p value
txt4 = uicontrol(fig_h, 'Style','text',...
    'Position',[pos(3)*margins pos(4)*0.06 pos(3)*0.15 pos(4)*0.1],...
    'String','p ANOVA', ...
    'HorizontalAlignment', 'left');
editPanova_h = uicontrol(fig_h, 'Style','edit',...
    'Position',[pos(3)*(margins+0.06) pos(4)*0.135 pos(3)*0.04 pos(4)*0.03], ...
    'String', ' < 0.01 ');

%------------------------- recovery probability
txt5 = uicontrol(fig_h, 'Style','text',...
    'Position',[pos(3)*margins pos(4)*0.01 pos(3)*0.15 pos(4)*0.1],...
    'String','p recov ', ...
    'HorizontalAlignment', 'left');

editPrecov_h = uicontrol(fig_h, 'Style','edit',...
    'Position',[pos(3)*(margins+0.06) pos(4)*0.085 pos(3)*0.04 pos(4)*0.03], ...
    'String', '<= 1 ');

%--------------------------- popup for Data Criteria - Highest dose, etc.
popC_h = uicontrol(fig_h, ....
    'Style',  'popupmenu',...
    'String', datacrit,...
    'Position', [pos(3)*margins pos(4)*0.7 pos(3)*0.15 pos(4)*0.1]);

%--------------------------- popup for Data Criteria - Electrode broken
popEcrit_h = uicontrol(fig_h, ....
    'Style',  'popupmenu',...
    'String', electrodecrit,...
    'Position', [pos(3)*margins pos(4)*0.65 pos(3)*0.15 pos(4)*0.1]);

%--------------------------- radiobuttons for monkeys
bg = uibuttongroup(...
    'Position', [0.03 0.82 0.15 0.15]);

% Create three radio buttons in the button group.
r1 = uicontrol(bg,'Style',...
    'radiobutton',...
    'String','both',...
    'Position',[1 5 100 15]);

r2 = uicontrol(bg,'Style','radiobutton',...
    'String','mango',...
    'Position',[1 35 100 15]);

r3 = uicontrol(bg,'Style','radiobutton',...
    'String','kaki',...
    'Position',[1 65 100 15]);

%---------------------------- plot button
uicontrol(fig_h,...
    'Style', 'pushbutton',...
    'Position', [pos(3)*0.4 pos(4)*margins pos(3)*0.08 pos(4)*0.05],...
    'String',   'plot', ...
    'Callback', @createPlot);

%---------------------------- plot and save figure 1
uicontrol(fig_h,...
    'Style', 'pushbutton',...
    'Position', [pos(3)*0.5 pos(4)*margins pos(3)*0.05 pos(4)*0.05],...
    'String',   'plot Fig1', ...
    'Callback', @plotFigure1);

%---------------------------- plot and save figure 2
uicontrol(fig_h,...
    'Style', 'pushbutton',...
    'Position', [pos(3)*0.56 pos(4)*margins pos(3)*0.05 pos(4)*0.05],...
    'String',   'plot Fig2', ...
    'Callback', @plotFigure2);

%---------------------------- plot and save figure 3
uicontrol(fig_h,...
    'Style', 'pushbutton',...
    'Position', [pos(3)*0.62 pos(4)*margins pos(3)*0.05 pos(4)*0.05],...
    'String',   'plot Fig3', ...
    'Callback', @plotFigure3);

%---------------------------- plot and save figure 5 
uicontrol(fig_h,...
    'Style', 'pushbutton',...
    'Position', [pos(3)*0.68 pos(4)*margins pos(3)*0.05 pos(4)*0.05],...
    'String',   'plot Fig5', ...
    'Callback', @plotFigure5);

%----------------------- close all button
uicontrol(fig_h, ...
    'Style', 'pushbutton', ...
    'Position', [pos(3)*0.75 pos(4)*margins pos(3)*0.08 pos(4)*0.05], ...
    'String', 'close all' ,...
    'Callback', @CloseAll);

%------------------- export figure button
uicontrol(fig_h,...
    'Style', 'pushbutton',...
    'Position', [pos(3)*0.85 pos(4)*margins pos(3)*0.08 pos(4)*0.05],...
    'String',   'export', ...
    'Callback', @exportPlot);

%------------------------- export barplot
uicontrol(fig_h,...
    'Style', 'pushbutton',...
    'Position', [pos(3)*0.95 pos(4)*margins pos(3)*0.04 pos(4)*0.05],...
    'String',   'bar xy', ...
    'Callback', @exportBar);


%------------------------ popup for stimuli condition
pop_Yspec = uicontrol(fig_h,...
    'Style', 'popupmenu',...
    'String', stimulicond,...
    'Position',[pos(3)*0.2 pos(4)*0.8 pos(3)*0.08 pos(4)*0.05]);

pop_Xspec = uicontrol(fig_h,...
    'Style','popupmenu',...
    'String', stimulicond,...
    'Position',[pos(3)*0.2 pos(4)*0.6 pos(3)*0.08 pos(4)*0.05]);


%------------------------ popup for eye dominance
pop_Xeye = uicontrol(fig_h,...
    'Style', 'popupmenu',...
    'String', eyeopts,...
    'Position',[pos(3)*0.2 pos(4)*0.75 pos(3)*0.08 pos(4)*0.05]);

pop_Yeye = uicontrol(fig_h,...
    'Style','popupmenu',...
    'String', eyeopts,...
    'Position',[pos(3)*0.2 pos(4)*0.55 pos(3)*0.08 pos(4)*0.05]);

%------------------------ figures to plot when data point click
% when you change the Strings here, you have to adapt them in DataPressed.m

fig2plot_check(1) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Tuning Curve', ...
    'Position', [pos(3)*0.2 pos(4)*0.10 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);

fig2plot_check(2) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Regression', ...
    'Position', [pos(3)*0.2 pos(4)*0.25 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);

fig2plot_check(3) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Recovery', ...
    'Position', [pos(3)*0.3 pos(4)*0.25 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);

fig2plot_check(4) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'LFP Gui', ...
    'Position', [pos(3)*0.2 pos(4)*0.20 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);

fig2plot_check(5) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Raster', ...
    'Position', [pos(3)*0.3 pos(4)*0.20 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);

fig2plot_check(6) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Wave Form', ...
    'Position', [pos(3)*0.2 pos(4)*0.15 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);

fig2plot_check(7) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Spike Density', ...
    'Position', [pos(3)*0.3 pos(4)*0.15 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);


fig2plot_check(8) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Phase Select.', ...
    'Position', [pos(3)*0.3 pos(4)*0.10 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);

fig2plot_check(9) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'smooth PSTH', ...
    'Position', [pos(3)*0.3 pos(4)*0.30 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);


fig2plot_check(10) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Variability', ...
    'Position', [pos(3)*0.2 pos(4)*0.30 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);


fig2plot_check(11) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Dummy', ...
    'Position', [pos(3)*0.3 pos(4)*0.35 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);


fig2plot_check(12) = uicontrol(fig_h, ...
    'Style', 'checkbox', ...
    'String', 'Dummy', ...
    'Position', [pos(3)*0.2 pos(4)*0.35 pos(3)*0.08 pos(4)*0.05], ...
    'Callback', @UpdateFigures);


%% publication plots
    function plotFigure1(~,~,~)
        % figure 1 - all plots that are necessary to generate Fig. 1
        % 5HT effect on spiking activity
        
        % prepare folders
        if exist(fullfile(cd, 'raw_figs'), 'dir') ~= 7
            mkdir('raw_figs');
            mkdir(fullfile('raw_figs\Orientation'));
            mkdir(fullfile('raw_figs\SpatialF'));
            mkdir(fullfile('raw_figs\Contrast'));
            mkdir(fullfile('raw_figs\Size'));
        end
        
        % plots showing mean firing rate
        addCross.Value = 0;             addUnity.Value = 1;
        addHistograms.Value = 0;        eqAxes.Value = 1;
                
        popX_h.Value= find(strcmp(axisopts, 'mean spike rate base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'mean spike rate drug'))+1;
               
        pop_Yspec.Value = 2; pop_Xspec.Value = 2;  % Orientation
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\MeanFiringRate.fig'); delete(h)
                
        pop_Yspec.Value = 3; pop_Xspec.Value = 3; % Spatial Frequency
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\MeanFiringRate.fig'); delete(h)
                
        pop_Yspec.Value = 4; pop_Xspec.Value = 4; % Contrast
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\MeanFiringRate.fig'); delete(h)
              
        pop_Yspec.Value = 5; pop_Xspec.Value = 5; % Size
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\MeanFiringRate.fig'); delete(h)
        
        
        % plots showing regression plots
        addCross.Value = 1;        addUnity.Value = 0;
        addHistograms.Value = 1;   eqAxes.Value = 0;
        editr2_h.String = '>0.7'; % the regression fit criteria
        
        popX_h.Value= find(strcmp(axisopts, 'gain change'))+1;
        popY_h.Value= find(strcmp(axisopts, 'additive change (rel)'))+1;
        
        pop_Yspec.Value = 2; pop_Xspec.Value = 2; 
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\Regression.fig'); delete(h)
        
        pop_Yspec.Value = 3; pop_Xspec.Value = 3;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\Regression.fig'); delete(h)
        
        pop_Yspec.Value = 4; pop_Xspec.Value = 4;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\Regression.fig'); delete(h)
        
        pop_Yspec.Value = 5; pop_Xspec.Value = 5;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\Regression.fig'); delete(h)
        
        % example units
        % tuning curves
        h = openfig(exinfo([exinfo.idi] == 91 & [exinfo.ocul]==-1).fig_tc);
        savefig(h, 'raw_figs\Orientation\Tc5HT.fig'); delete(h);
        h = openfig(exinfo([exinfo.idi] == 412 & [exinfo.ocul]==1).fig_tc);
        savefig(h, 'raw_figs\Orientation\TcNaCl.fig'); delete(h);
        
        h = openfig(exinfo([exinfo.idi] == 124 & [exinfo.ocul]==-1).fig_tc);
        savefig(h, 'raw_figs\SpatialF\Tc5HT.fig'); delete(h);
        h = openfig(exinfo([exinfo.idi] == 627 & [exinfo.ocul]==-1).fig_tc);
        savefig(h, 'raw_figs\SpatialF\TcNaCl.fig'); delete(h);
        
        h = openfig(exinfo([exinfo.idi] == 510 & [exinfo.ocul]==1).fig_tc);
        savefig(h, 'raw_figs\Contrast\Tc5HT.fig'); delete(h);
        h = openfig(exinfo([exinfo.idi] == 492 & [exinfo.ocul]==-1).fig_tc);
        savefig(h, 'raw_figs\Contrast\TcNaCl.fig'); delete(h);
        
        h = openfig(exinfo([exinfo.idi] == 581 & [exinfo.ocul]==1).fig_tc);
        savefig(h, 'raw_figs\Size\Tc5HT.fig'); delete(h);
        h = openfig(exinfo([exinfo.idi] == 573 & [exinfo.ocul]==1).fig_tc);
        savefig(h, 'raw_figs\Size\TcNaCl.fig'); delete(h);
        
        %regression plots
        h = openfig(exinfo([exinfo.idi] == 91 & [exinfo.ocul]==-1).fig_regl);
        savefig(h, 'raw_figs\Orientation\Reg5HT.fig'); delete(h);
        h = openfig(exinfo([exinfo.idi] == 412 & [exinfo.ocul]==1).fig_regl);
        savefig(h, 'raw_figs\Orientation\RegNaCl.fig'); delete(h);
        
        h = openfig(exinfo([exinfo.idi] == 124 & [exinfo.ocul]==-1).fig_regl);
        savefig(h, 'raw_figs\SpatialF\Reg5HT.fig'); delete(h);
        h = openfig(exinfo([exinfo.idi] == 627 & [exinfo.ocul]==-1).fig_regl);
        savefig(h, 'raw_figs\SpatialF\RegNaCl.fig'); delete(h);
        
        h = openfig(exinfo([exinfo.idi] == 510 & [exinfo.ocul]==1).fig_regl);
        savefig(h, 'raw_figs\Contrast\Reg5HT.fig'); delete(h);
        h = openfig(exinfo([exinfo.idi] == 492 & [exinfo.ocul]==-1).fig_regl);
        savefig(h, 'raw_figs\Contrast\RegNaCl.fig'); delete(h);
        
        h = openfig(exinfo([exinfo.idi] == 581 & [exinfo.ocul]==1).fig_regl);
        savefig(h, 'raw_figs\Size\Reg5HT.fig'); delete(h);
        h = openfig(exinfo([exinfo.idi] == 573 & [exinfo.ocul]==1).fig_regl);
        savefig(h, 'raw_figs\Size\RegNaCl.fig'); delete(h);
        %%
        
    end

    function plotFigure2(~,~,~)
        % figure 2 - all plots that are necessary to generate Fig. 2
        % 5HT effect on tuning curves
    
        % prepare folders
        if exist('raw_figs', 'dir') ~= 7
            mkdir('raw_figs');
            mkdir(fullfile('raw_figs\Orientation'));
            mkdir(fullfile('raw_figs\SpatialF'));
            mkdir(fullfile('raw_figs\Contrast'));
            mkdir(fullfile('raw_figs\Size'));
        end
        
        % plots showing parameters of the best fitting tuning function        
        addCross.Value = 0;         addUnity.Value = 1;
        addHistograms.Value = 0;    eqAxes.Value = 1;        
        editr2tc_h.String = '>0.7'; % criteria of the tuning curve fit
        
        % ORIENTATION
        pop_Yspec.Value = 2; pop_Xspec.Value = 2;
        
        popX_h.Value= find(strcmp(axisopts, 'gauss fit mu base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss fit mu drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\PreferredOR.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss fit sig base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss fit sig drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\BandWidth.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss fit a base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss fit a drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\Amplitude.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss fit b base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss fit b drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\Offset.fig'); delete(h)
        
        
        % SPATIAL FREQUENCY
        pop_Yspec.Value = 3; pop_Xspec.Value = 3;
        
        popX_h.Value= find(strcmp(axisopts, 'gauss fit mu base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss fit mu drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\PreferredSF.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss fit sig base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss fit sig drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\BandWidth.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss fit a base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss fit a drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\Amplitude.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss fit b base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss fit b drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\Offset.fig'); delete(h)
        
        % CONTRAST
        pop_Yspec.Value = 4; pop_Xspec.Value = 4;
        
        popX_h.Value= find(strcmp(axisopts, 'c50 base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'c50 drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\C50.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'rmax base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'rmax drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\Rmax.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'co fit n base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'co fit n drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\ExponentN.fig'); delete(h)
        
        
        popX_h.Value= find(strcmp(axisopts, 'co fit m base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'co fit m drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\SpontActivity.fig'); delete(h)
        
        % SIZE
        pop_Yspec.Value = 5; pop_Xspec.Value = 5;
        
        popX_h.Value= find(strcmp(axisopts, 'SI base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'SI drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\SuppressionIndex.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'pref size base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'pref size drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\PreferredSize.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss ratio fit width center'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss ratio fit width center drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\WdtCenter.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss ratio fit width surround'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss ratio fit width surround drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\WdtSurround.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss ratio fit gain surround'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss ratio fit gain surround drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\GainSurround.fig'); delete(h)
        
        popX_h.Value= find(strcmp(axisopts, 'gauss ratio fit gain center'))+1;
        popY_h.Value= find(strcmp(axisopts, 'gauss ratio fit gain center drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\GainCenter.fig'); delete(h)
        
    end

    function plotFigure3(~,~,~)
        %figure 3 - all plots necessary to generate Fig. 3
        % 5HT effect on response dynamics
        
        pop_Yspec.Value = 6; pop_Xspec.Value = 6;
        
        % prepare folders
        if exist('raw_figs', 'dir') ~= 7
            mkdir('raw_figs');
            mkdir(fullfile('raw_figs\RC_population'));
        end
        
        % plots showing mean firing rates
        addCross.Value = 0; addHistograms.Value = 0;
        addUnity.Value = 1; eqAxes.Value = 1; addRegress.Value = 0;

        editr2_h.String = '>-inf';
        editr2tc_h.String = '>-inf';
        editlatency_h.String = '>-inf';
        
        
        popX_h.Value= find(strcmp(axisopts, 'mean spike rate base'))+1;
        popY_h.Value= find(strcmp(axisopts, 'mean spike rate drug'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\RC_population\MeanFiringRate.fig'); delete(h)
        
        
        % plots showing response latencies
        editlatency_h.String = '>0';
        
        popX_h.Value= find(strcmp(axisopts, 'latency base corrected'))+1;
        popY_h.Value= find(strcmp(axisopts, 'latency drug corrected'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\RC_population\Latency.fig'); delete(h)
        
        
        % plots showing latency change vs relative firing rate
        addCross.Value = 1; addUnity.Value = 0; eqAxes.Value = 0;
        addRegress.Value = 1;
        editlatency_h.String = '>0';

        popX_h.Value= find(strcmp(axisopts, 'nonparam area ratio'))+1;
        popY_h.Value= find(strcmp(axisopts, 'latency diff corrected'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\RC_population\Corr_Lat_RelRate.fig'); delete(h)
        
        
        % plot showing type-II regression results
        addCross.Value = 1;  addHistograms.Value = 1;
        addUnity.Value = 0;  eqAxes.Value = 0;
        addRegress.Value = 0;
        
        editlatency_h.String = '>-inf';
        editr2_h.String = '>0.7';
        
        popX_h.Value= find(strcmp(axisopts, 'gain change'))+1;
        popY_h.Value= find(strcmp(axisopts, 'additive change (rel)'))+1;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\RC_population\Regression.fig'); delete(h)
        
    end

    function plotFigure5(~,~,~)
        %figure 6 - all plots necessary to generate Fig. 6
        
        % determine which configuration is supposed to be shown, i.e.
        % whether to use the analysis based on the full experiments or the
        % analysis based on only the second part of the experiment
        suffix = ' 2nd half';  % '';
        
        
        
        % for the noise correlation choose whether you want to see the
        % simple or the spike count corrected metric 
        suffix2 = ' corrected'; 
        
        % prepare folders
        if exist('raw_figs', 'dir') ~= 7
            mkdir('raw_figs');
            mkdir(fullfile('raw_figs\Orientation'));
            mkdir(fullfile('raw_figs\SpatialF'));
            mkdir(fullfile('raw_figs\Contrast'));
            mkdir(fullfile('raw_figs\Size'));
        end
        
        % plots showing variability and co-variability
        addCross.Value = 1; addHistograms.Value = 0;
        addUnity.Value = 1; eqAxes.Value = 1; addRegress.Value = 0;

        editr2_h.String = '>-inf';
        editr2tc_h.String = '>-inf';
        editlatency_h.String = '>-inf';
        
        % fano factor
        popX_h.Value= find(strcmp(axisopts, ['fano factor' suffix ' base']))+1;
        popY_h.Value= find(strcmp(axisopts, ['fano factor' suffix ' drug']))+1;

        pop_Yspec.Value = 2; pop_Xspec.Value = 2;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\FanoFactor.fig'); delete(h)
        
        pop_Yspec.Value = 3; pop_Xspec.Value = 3;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\FanoFactor.fig'); delete(h)
        
        pop_Yspec.Value = 4; pop_Xspec.Value = 4;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\FanoFactor.fig'); delete(h)
        
        pop_Yspec.Value = 5; pop_Xspec.Value = 5;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\FanoFactor.fig'); delete(h)
        
        
        % fano factor vs relative rate
        popX_h.Value= find(strcmp(axisopts, ['nonparam area ratio' suffix]))+1;
        popY_h.Value= find(strcmp(axisopts, ['fano factor' suffix ' diff']))+1;
        addUnity.Value = 0;  eqAxes.Value = 0;

        pop_Yspec.Value = 2; pop_Xspec.Value = 2;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\FanoFactor_vs_RelRate.fig'); delete(h)
        
        pop_Yspec.Value = 3; pop_Xspec.Value = 3;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\FanoFactor_vs_RelRate.fig'); delete(h)
        
        pop_Yspec.Value = 4; pop_Xspec.Value = 4;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\FanoFactor_vs_RelRate.fig'); delete(h)
        
        pop_Yspec.Value = 5; pop_Xspec.Value = 5;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\FanoFactor_vs_RelRate.fig'); delete(h)
        
        
        % noise correlation
        popX_h.Value= find(strcmp(axisopts, ['noise correlation' suffix suffix2]))+1;
        popY_h.Value= find(strcmp(axisopts, ['noise correlation drug' suffix suffix2]))+1;
        addUnity.Value = 1;  eqAxes.Value = 1;

        pop_Yspec.Value = 2; pop_Xspec.Value = 2;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\NoiseCorr.fig'); delete(h)
        
        pop_Yspec.Value = 3; pop_Xspec.Value = 3;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\NoiseCorr.fig'); delete(h)
        
        pop_Yspec.Value = 4; pop_Xspec.Value = 4;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\NoiseCorr.fig'); delete(h)
        
        pop_Yspec.Value = 5; pop_Xspec.Value = 5;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\NoiseCorr.fig'); delete(h)
       
        
        % noise correlation vs relative rate
        popX_h.Value= find(strcmp(axisopts, ['nonparam area ratio' suffix ]))+1;
        popY_h.Value= find(strcmp(axisopts, ['noise correlation diff' suffix suffix2]))+1;
        addUnity.Value = 0;  eqAxes.Value = 0;

        pop_Yspec.Value = 2; pop_Xspec.Value = 2;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Orientation\NoiseCorr_vs_RelRate.fig'); delete(h)
        
        pop_Yspec.Value = 3; pop_Xspec.Value = 3;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\SpatialF\NoiseCorr_vs_RelRate.fig'); delete(h)
        
        pop_Yspec.Value = 4; pop_Xspec.Value = 4;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Contrast\NoiseCorr_vs_RelRate.fig'); delete(h)
        
        pop_Yspec.Value = 5; pop_Xspec.Value = 5;
        createPlot(0, 0, 0); h = exportPlot(0, 0);
        savefig(h, 'raw_figs\Size\NoiseCorr_vs_RelRate.fig'); delete(h)
    end



%% nested functions
    
    function exportBar(~, ~)
        % creates a bar plot that compares the x and y values for each drug
        % condition.
        
        figure;
        
        if any(dat.is5HT) && any(~dat.is5HT)
            % considering a bar plot for each drug condition 
            n_5HT = sum(dat.is5HT);
            n_nacl = sum(~dat.is5HT);
            
            bar([1 2 4 5], [mean(dat.y(dat.is5HT)), mean(dat.x(dat.is5HT)),...
                mean(dat.y(~dat.is5HT)) mean(dat.x(~dat.is5HT))]); hold on;
            errorbar([1 2 4 5], [mean(dat.y(dat.is5HT)), mean(dat.x(dat.is5HT)),...
                mean(dat.y(~dat.is5HT)), mean(dat.x(~dat.is5HT))], ...
                [std(dat.y(dat.is5HT))/sqrt(n_5HT), std(dat.x(dat.is5HT))/sqrt(n_5HT),...
                std(dat.y(~dat.is5HT))/sqrt(n_nacl), std(dat.x(~dat.is5HT))/sqrt(n_nacl)], '.');
            
            set(gca, 'XLim', [0.5 5.5], 'XTick', [1 2 4 5], ...
                'XTickLabel', {'5HT' , '5HT base', ...
                'NaCl', 'NaCl base'});
            ylabel(dat.xlab);
            
            p_wil_5HT = ranksum(dat.x(dat.is5HT), dat.y(dat.is5HT));
            [~, p_ttest_5HT] = ttest(dat.x(dat.is5HT), dat.y(dat.is5HT));
            
            p_wil_nacl = ranksum(dat.x(~dat.is5HT), dat.y(~dat.is5HT));
            [~, p_ttest_nacl] = ttest(dat.x(~dat.is5HT), dat.y(~dat.is5HT));
            
            p_wil_vs = ranksum(dat.y(dat.is5HT), dat.y(~dat.is5HT));
            [~, p_ttest_vs] = ttest2(dat.y(dat.is5HT), dat.y(~dat.is5HT));
            
            title(sprintf(['5HT (n=%1.0f): ranksum p=%1.3f \t ttest p=%1.3f \n' ...
                'NaCl (n=%1.0f): ranksum p=%1.3f \t ttest p=%1.3f \n' ...
                '5HT vs NaCl (scatter y-axis): ranksum p=%1.3f \t ttest2 p=%1.3f \n \n' ...
                '5HT mean %1.2f, mean control %1.2f \t \t NaCl mean %1.2f, mean control %1.2f'], ...
                n_5HT, p_wil_5HT, p_ttest_5HT, ...
                n_nacl, p_wil_nacl, p_ttest_nacl, ...
                p_wil_vs, p_ttest_vs, ...
                mean(dat.y(dat.is5HT)), mean(dat.x(dat.is5HT)),...
                mean(dat.y(~dat.is5HT)), mean(dat.x(~dat.is5HT))));
        else
            % assuming all data 
            n = length(dat.x);
            nsqr = sqrt(n);
            bar(1, mean(dat.y), 'FaceColor', 'r', 'EdgeColor', 'r', 'LineWidth', 2); hold on;
            bar(2, mean(dat.x), 'FaceColor', 'w', 'EdgeColor', 'r'); hold on;
            errorbar([1 2], [mean(dat.y), mean(dat.x)], ...
                [std(dat.y)/nsqr, std(dat.x)/nsqr], 'k.', 'LineWidth', 2); hold on;
            plot([0.5 2.5], [0 0], 'k');
            set(gca, 'XLim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel', {'Serotonin', 'Baseline'});
            ylabel(dat.xlab);
            
            p_wil = ranksum(dat.x, dat.y);
            [~, p_ttest] = ttest(dat.x, dat.y);
            
            if any(dat.is5HT)
                title(sprintf(['5HT (n=%1.0f) mu _{5HT}=%1.3f, med_{5HT}=%1.3f \n'...
                    'mu _{Base}=%1.3f med_{Base}=%1.3f \n'...
                    'ranksum p=%1.3f \t ttest p=%1.3f'], ...
                    n, mean(dat.y), median(dat.y), mean(dat.x), median(dat.x), p_wil, p_ttest));
            else
                title(sprintf('NaCl (n=%1.0f) ranksum p=%1.3f \t ttest p=%1.3f', n, p_wil, p_ttest));
            end
            
        end
        set(gcf, 'UserData', dat);
        hold off;
        
    end

    function h2 = exportPlot(~, ~)
        % copys the plot with all histograms and annotations into a figure
        % outside the gui.
        % it also addts the statistics to the dat structure and attaches
        % this as user data to the figure
        
        h2 = figure();
        ax = findobj(fig_h, 'type', 'axes');
        if length(ax) == 1
            ax_new = copyobj(ax, h2, 'legacy');
            ax_new.Position = [ 0.15 0.15 0.7 0.7];
        elseif length(ax) == 3
            ax_new = copyobj(ax, h2, 'legacy');
            
            ax_new(1).Position = [0.15 0.1 0.45 0.55];
            ax_new(2).Position = [0.7 0.1 0.15 0.55];
            ax_new(3).Position = [0.15 0.7 0.45 0.15];
            h2.Position = [643 328 737 537];
        else
            warning('there are too many axes to handle');
        end
        
        set(h2, 'UserData', dat);
        setStats( h2, getInclusionCrit() );
        
    end

% -------------------------------------------------------------------------

    function createPlot(~, eventdata, ~)
        % plot the analysis results according to the options entered via the GUI
        
        if (get(popY_h, 'Value') == 1) || (get(popX_h, 'Value') == 1) % no axis optiosn were chosen
            allax = findobj(gcf,'type','axes'); 
            delete(allax);
        else
            allax = findobj(gcf,'type','axes');
            delete(allax);
            
            fctX = axisopts{get(popX_h, 'Value')-1};
            fctY = axisopts{get(popY_h, 'Value')-1};
            spec.stimx = stimulicond{get(pop_Xspec, 'Value')};
            spec.stimy = stimulicond{get(pop_Yspec, 'Value')};
            
            spec.eyex = eyeopts{get(pop_Xeye, 'Value')};
            spec.eyey = eyeopts{get(pop_Yeye, 'Value')};
            
            
            % special case for mitchel plot 
            % instead of the scatter plot, the average variance for binned
            % mean responses is shown
            if contains(fctY, 'mitchel') && ~contains(spec.stimy, 'all')
                UpdateInclusion(editlatency_h, eventdata);
                
                % only 5HT experiments and non RC experiments
                incl_i =  intersect( incl_i, find(~[exinfo.isRC] & [exinfo.is5HT])); 
                createPlotHelper(get(addHistograms, 'Value'), spec);
                
                %dat.exinfo = exinfo(incl_i);
                set(gcf, 'UserData', dat);
                
            else
                
                UpdateInclusion(editlatency_h, eventdata);
                dat = createUnitPlot(exinfo(incl_i), fctX, fctY, spec, ...
                    suFigsOnDataClick, get(addHistograms, 'Value'));
                
                %dat.exinfo = exinfo(incl_i);
                set(gcf, 'UserData', dat);
                
            end
            
            UpdateAxes(fctX, fctY);
            
        end
        
        if eqAxes.Value;                eqax;        end
        if addUnity.Value;             unity;        end
        if addCross.Value;            crossl;        end
        if addRegress.Value;            regl;        end
        
    end

%-------------------------------------------------------------------------
    function UpdateInclusion(~, ~)
        % This functions tries to exclude all data that are note relevant
        % for the plot. This includes general inclusion criteria such as
        % the number of repeats for each stimulus, the electrode resistance
        % or the minimal firing rate.
        
        
        incl_i = 1:length(exinfo);  % in the beginning assume all experiment are included
        
        % stimulus conditions
        stimy = pop_Yspec.String{pop_Yspec.Value};
        stimx = pop_Xspec.String{pop_Xspec.Value};
        
        % which stimuli condition should be shown
        UpdateInclusionHelper( ...
            strcmp(stimy, stimx) & ~strcmp(stimy, 'RC') & ~strcmp(stimy, 'all stimuli cond'), ...
            ['cellfun(@(x) strcmp(x, ''' stimx '''), {exinfo.param1})']);
        
        UpdateInclusionHelper( strcmp(stimy, stimx) & ...
            ~strcmp(stimy, 'RC')  & strcmp(stimy, 'all stimuli cond'), ...
            '~[exinfo.isRC]');

        
        % the electrode resistance
%         UpdateInclusionHelper(1,...
%             'cellfun(@(x) x<150, {exinfo.resistance})');
        UpdateInclusionHelper(1,...
            'cellfun(@(x) isnan(x) || x<150, {exinfo.resistance})');
        
        
        %%% global inclusion criteria
        % the minimum firing rate for the preferred stimulus
        UpdateInclusionHelper(~strcmp(stimy, 'RC'),...
            ['cellfun(@max, {exinfo.ratemn}) >=' num2str(minspk)]);
        UpdateInclusionHelper(strcmp(stimy, 'RC'),...
            ['cellfun(@max, {exinfo.ratemn}) >=' num2str(minspk/100)]);
        
        % number of repetitions greater than the threshold set in the first
        % code lines
        UpdateInclusionHelper(1,...
            ['cellfun(@min, {exinfo.nrep}) >=' num2str(nrep)]);
        UpdateInclusionHelper(1,...
            ['cellfun(@min, {exinfo.nrep_drug}) >=' num2str(nrep)]);
        
        
        %%% GUI set inclusion criteria
        % which eye criteria should be set
        UpdateInclusionHelper(strcmp(pop_Xeye.String(pop_Xeye.Value), 'dominant eye'), '[exinfo.isdominant]');
        UpdateInclusionHelper(strcmp(pop_Xeye.String(pop_Xeye.Value), 'non-dominant eye'), '~[exinfo.isdominant]');
        
        
        % which monkey data should be shown
        UpdateInclusionHelper(r2.Value, '[exinfo.ismango]');
        UpdateInclusionHelper(r3.Value, '~[exinfo.ismango]');
        
        
        % what is the regression fit threshold
        UpdateInclusionHelper(1, ['[exinfo.r2reg] ' get(editr2_h, 'String')]);
        
        
        % what is the latency threshold
        datalatency = evalMU('latency base', 'latency drug', exinfo);
        if strcmp(stimulicond(get(pop_Xspec, 'Value')), 'RC')
            cond =  ['datalatency.x'  get(editlatency_h, 'String') ' & ' ...
                'datalatency.y'  get(editlatency_h, 'String') ];
            UpdateInclusionHelper(true, cond)
        end
        
        % what is the tuning fit threshold
        UpdateInclusionHelper(1,...
            ['[exinfo.gaussr2] ' get(editr2tc_h, 'String') ...
            ' & [exinfo.gaussr2_drug] ' get(editr2tc_h, 'String')])
        
        % what is the anova alpha level
        UpdateInclusionHelper(1,...
            ['[exinfo.p_anova] ' get(editPanova_h, 'String') ...
            ' & [exinfo.p_anova_drug] ' get(editPanova_h, 'String')])
        
        % what is the isolation quality criteria
        UpdateInclusionHelper(1,...
            ['[exinfo.spkqual_base] ' get(editSpkSortQualtiy_h, 'String') ...
            ' & [exinfo.spkqual_drug] ' get(editSpkSortQualtiy_h, 'String')])

        
        % which one of the experiments in one unit recording should be shown
        idx1 = getCritIdx(exinfo(incl_i), popC_h.String{popC_h.Value});
        incl_i = incl_i(idx1);
        
        % which electrode criteria is considered
        idx2 = getElectrodeCrit(exinfo(incl_i), popEcrit_h.String{popEcrit_h.Value});
        incl_i = incl_i(idx2);
        
        
        %in case you want to exclude a single data point
        incl_i( ismember(incl_i, find([ exinfo.id ] == 30) )) =[]; % mango#30 - bad spike sorting
        incl_i( ismember(incl_i, find([ exinfo.id ] == 65) )) = []; % mango#65 - probably in V2


        % for variability metrices concatenated experiments should be
        % excluded because they show some adaptation in the beginning of
        % each individual experiment that artificially increases the
        % variability metrices
        fctX = axisopts{get(popX_h, 'Value')-1};
        fctY = axisopts{get(popY_h, 'Value')-1};
        if ( contains(fctX, 'fano factor') || contains(fctY, 'fano factor') || ...
            contains(fctX, 'noise correlation') || contains(fctY, 'noise correlation') ) ...
            && (contains(fctX, '2nd half') || contains(fctY, '2nd half'))
            incl_i = intersect(incl_i, ...
            find( ~contains({exinfo.fname}, 'all.grating') & ...
                ~contains({exinfo.fname_drug}, 'all.grating') ));
        end
        
        
        % what is the threshold for the recovery p-value
        UpdateInclusionHelper(1, ['[exinfo.ret2base] ' get(editPrecov_h, 'String')]);
        
    end

    function UpdateInclusionHelper(eval_flag, cond)
        if eval_flag
            eval(['incl_i = intersect( incl_i, find(' cond '));']);
            fprintf([cond '\n']);
        end
    end
%--------------------------------------------------------------------------
    function UpdateAxes(fctX, fctY)
        % update the axes according to the plotted values, especially to
        % transform the axis scale to logarithmic scale and to predefine
        % the limits of the axes
        
        
        axrg = [0.04 64];
        allax = findobj(gcf, 'Type', 'axes');
        
        % X AXES
        if (contains(fctX, 'gain') && ~contains(fctX, 'fit')) || ...
                contains(fctX, 'nonparam')
            
            % turn to log scaled axes
            if length(allax) > 1
                plotHist(dat.x, dat.is5HT, allax(3), 'log');
                set(allax(3), 'xlim', axrg, 'XScale','log', 'XTick', [0.25 1 4]);
            end
            set(allax(1), 'xlim', axrg, 'XScale','log', 'XTick', [0.25 1 4]);
        end
        
        % Y AXES
        if (contains(fctY, 'gain') && ~contains(fctY, 'fit'))...
                || contains(fctY, 'nonparam') 
            
            % turn to log scaled axes
            if length(allax) >1
                plotHist(dat.y, dat.is5HT, allax(2), 'log');
                set(allax(2), 'XScale','log', 'XTick', [0.25 1 4], 'xlim', axrg);
            end
            
            set(allax(1), 'YScale','log', 'YTick', [0.25 1 4], 'ylim', axrg);
        end
    end

%--------------------------------------------------------------------------
    function UpdateFigures(~, eventdata, ~)
        % add the clicked checkbox to the figures that should pop up upon
        % clicking of a data point in the plot
        
        suFigsOnDataClick = vertcat({fig2plot_check.String});
        suFigsOnDataClick = suFigsOnDataClick([fig2plot_check.Value]==1);
        
        % update the plot because each data point has to have the new
        % figure names in its ButtonDownFcn
        createPlot(0, eventdata, 0);
    end

%--------------------------------------------------------------------------
    function s = getInclusionCrit()
        % concatenate all inclusion criteria and their values set by the
        % GUI to a long string. the string is attached to the figure and
        % makes it easier to later comprehend how the plot was generated.
        
        % all inclusion criteria as string
        s = sprintf([popC_h.String{popC_h.Value} '  ' popEcrit_h.String{popEcrit_h.Value} '\n' ...
            'min # trial repetitions ' num2str(nrep) '\n'...
            'min # spikes in preferred condition ' num2str(minspk) '\n' ...
            'anova p ' get(editPanova_h, 'String') '\n' ...
            'overlap between baseline and control bootstrap p ' get(editPrecov_h, 'String') '\n' ...
            'regression fit r2 ' get(editr2_h, 'String') '\n' ...
            'tuning curve fit r2 ' get(editr2tc_h, 'String') '\n' ...
            'latency estimation ' get(editlatency_h, 'String') '\n' ...
            'spike isolation quality ' get(editSpkSortQualtiy_h, 'String')]);
        
    end

end

%--------------------------------------------------------------------------
function CloseAll(~, ~, ~)
    % close all open figures
    close all;
end


function my_closereq(~,~)
% Close request function
% to display a question dialog box upon closing of a figure

selection = questdlg('Close This Figure?',...
    'Close Request Function',...
    'Yes','No','Yes');

switch selection
    case 'Yes'
        delete(gcf)
    case 'No'
        return
end
end



function idx = getElectrodeCrit(exinfo, crit)
% returns the index for exinfo that fullfiles the electrode criteria crit


switch crit
    case 'all'
        idx = 1:length(exinfo);
    case 'unbroken and included broken ones'
        idx = find(~[exinfo.electrodebroken_excl]);
    case 'unbroken only'
        idx = find(~[exinfo.electrodebroken]);
    case 'included broken ones'
        idx = find([exinfo.electrodebroken] & ~[exinfo.electrodebroken_excl]);
    case 'broken only'
        idx = find([exinfo.electrodebroken]);
    case 'neg volt'
        idx = find([exinfo.volt]<0);
end

end




