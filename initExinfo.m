function exinfo = initExinfo(varargin)
% exinfo = initExinfo()
%
% This function initializes the result structure exinfo with each row
% dedicated to the analysis of two experiments and their differences. More
% specifically, each row contains the results of the analysis of each
% individual experiment (i.e. spike rate, varaibility, tuning, etc.) and
% their comparison (i.e. regression fit, relative activity change, etc.).  
% 
% % This function also adds...
% - all the generic background information (mostly from the filename or ex file), 
%     e.g. dose, drug, spike waveform cluster, etc.
% - information about the electrode (hardcoded within the code),
%     e.g. was it broken, is it still considere for analysis
% - and specifies the folder and names of figures (also hardcoded)  
%     e.g. the figure name of the tuning curve, regression fit, etc.
% 
% The output (exinfo) is stored in a superordinate folder called Data as
% 'empty_exinfo.mat'. If the file already exists, it is overwritten.
%
%
%
% Optional arguments:
% 'fname', XX     - XX being the filename of the text file containing the
%                   pairs of experiments. If not specified, it is 
%                   'Z:\Corinna\filenames\SU110716_CL_all.txt'
% 'exinfo', XX    - XX being the result structure. Note that the algorithm
%                   still loops through the filenames in fname.
%
% 'figdir', XX    - XX being a folder name, destination for all matlab
%                   figures associated with this exinfo structure
%
% 
% The function further calls:
% - setReceptiveFieldSize   - determines the RF width by finding and
%                               analyzing the YPos.mat and XPos.mat files
% - addSortingValue         - retrieves the spike sorting quality
% - addNumInExp             - retrieves information about the preceeding
%                               drug experiments
%  
% 
% 
% 
% @CL 
% commented 19.05.2017 

%% add all subfolders to the path
addpath(genpath(pwd)); 
addpath(genpath('C:\Users\katsuhisa\Documents\code\analysis\integrated\interaction_project\GenAlyz_Code'));
addpath(genpath('C:\Users\katsuhisa\Documents\code\analysis\integrated\corinnas\HNetc'));

%% initiate variables
kk = 0;
idi = 1;
fname = 'Z:\Corinna\filenames\SU_CL_all.txt'; % txt file containing the experiment filenames

% mypath = 'D:\Users\kk\interaction_project\dataset\';
mypath = 'Z:\Katsuhisa\interaction_project\dataset\';
figdir.pre = fullfile(mypath, 'Figures\'); %folder destination for all figures
figdir.Data = fullfile(mypath, 'Data\');%folder destination for data


% result structure
exinfo = struct(...  
    'id', [], ... %<- unit ID. Mango is indicated by ID.0, Kaki by ID.5
    'idi', [], ...  %<- in case there are multiple ocular conditions, they all have the same idi. This is to identify them.
    'monkey', {}, ... %<- string indicating mango or kaki
    'fname', [], 'fname_drug', [], ... %<- filename of the baseline and drug ex file
    'date', [], 'date_drug', [], ... %<- date of the recording
    'ocul', [], ... %<- occularity
    'drugname', [], ... %<- 5HT or NaCl
    'dose', [], ... %<- applied iontophoresis current (nA), see getDose(ex)
    'dosernd', [], ... %<- binned dose, see getDose(ex)
    'volt', [],... %<- recorded voltage, see getVolt(ex)
    'resistance', [], ... %<- computed as volt/dose
    'gslope', [], ... %<- type-II regression slope
    'yoff', [], ... %<- type-II regression y offset
    'r2reg', [], ... %<- type-II regression explained variance
    'yoff_rel', [], ... %<- type-II regression y offset for normalized spike rates (note that the slope and fit do not change)
    'param1', {},... %<- first stimulus dimension tested (or, co, sz, ...)
    'param2', {}, ... %<- second stimulus dimension tested (me or co), important for RC orxco files
    'nonparam_ratio', [], ... %<- ratio between mean firing rate
    'lat', [], 'lat_drug', [], ... %<- latency estimate by maximum-likelihood method
    'lat2Hmax', [], 'lat2Hmax_drug', [], ... %<- latency estimates as time at half-maximum response
    'dur', [], 'dur_drug', [], ... %<- response duration, i.e. time between first and second half-maximum response
    'lat_1', [], 'lat_1_drug', [], ... %<- latency estimate by maximum-likelihood method
    'lat2Hmax_1', [], 'lat2Hmax_1_drug', [], ... %<- latency estimates as time at half-maximum response
    'dur_1', [], 'dur_1_drug', [], ... %<- response duration, i.e. time between first and second half-maximum response
    'lat_2', [], 'lat_2_drug', [], ... %<- latency estimate by maximum-likelihood method
    'lat2Hmax_2', [], 'lat2Hmax_2_drug', [], ... %<- latency estimates as time at half-maximum response
    'dur_2', [], 'dur_2_drug', [], ... %<- response duration, i.e. time between first and second half-maximum response
    'fitparam', [], 'fitparam_drug', [], ... %<- tuning fit parameters (e.g. orientation: mu, sigma, amplitude and offset)
    'fitparam_1', [], 'fitparam_1_drug', [], ... %<- tuning fit parameters (e.g. orientation: mu, sigma, amplitude and offset)
    'fitparam_2', [], 'fitparam_2_drug', [], ... %<- tuning fit parameters (e.g. orientation: mu, sigma, amplitude and offset)
    'ratemn' , [], 'ratemn_drug', [],... %<- mean firing rate for each stimulus condition
    'ratevars', [], 'ratevars_drug', [], ...%<- firing rate variance for each stimulus condition
    'ratesme', [], 'ratesme_drug', [],... %<- SEM of each firing rate
    'ratesd', [], 'ratesd_drug', [],... %<- SD of each firing rate
    'ratemn_1' , [], 'ratemn_1_drug', [],... %<- mean firing rate for each stimulus condition (pupil median split: small)
    'ratevars_1', [], 'ratevars_1_drug', [], ...%<- firing rate variance for each stimulus condition (pupil median split: small)
    'ratesme_1', [], 'ratesme_1_drug', [],... %<- SEM of each firing rate (pupil median split: small)
    'ratesd_1', [], 'ratesd_1_drug', [],... %<- SD of each firing rate (pupil median split: small)
    'ratemn_2' , [], 'ratemn_2_drug', [],... %<- mean firing rate for each stimulus condition (pupil median split: large)
    'ratevars_2', [], 'ratevars_2_drug', [], ...%<- firing rate variance for each stimulus condition (pupil median split: large)
    'ratesme_2', [], 'ratesme_2_drug', [],... %<- SEM of each firing rate (pupil median split: large)
    'ratesd_2', [], 'ratesd_2_drug', [],... %<- SD of each firing rate (pupil median split: large)
    'psglm',[], ... %<- glm data (b: beta weight for drug, pupil, their interaction; varexp: variance explained for stm only, + drug, + pupil, + interaction)
    'b_pf_drug', [],  ...%<- weight for drug to predict spike rate in preferred stimulus
    'b_pf_pupil', [], ... %<- weight for pupil to predict spike rate in preferred stimulus
    'b_pf_intr', [], ... %<- weight for interaction to predict spike rate in preferred stimulus
    'b_upf_drug', [],  ...%<- weight for drug to predict spike rate in unpreferred stimulus
    'b_upf_pupil', [], ... %<- weight for pupil to predict spike rate in unpreferred stimulus
    'b_upf_intr', [], ... %<- weight for interaction to predict spike rate in unpreferred stimulus
    'rawspkrates', [], 'rawspkrates_drug', [],... %<- ???
    'rate_resmpl', [], 'rate_resmpl_drug', [],... <- resampled firing rates
    'ratepar', [], 'ratepar_drug', [], ... %<- corresponding stimulus parameters
    'nrep', [], 'nrep_drug', [], ... %<- number of repeats per stimulus condition (only rewarded trials)
    'rawspkrates_1', [], 'rawspkrates_1_drug', [],... %<- ???
    'rate_resmpl_1', [], 'rate_resmpl_1_drug', [],... <- resampled firing rates
    'ratepar_1', [], 'ratepar_1_drug', [], ... %<- corresponding stimulus parameters
    'nrep_1', [], 'nrep_1_drug', [], ... %<- number of repeats per stimulus condition (only rewarded trials)
    'rawspkrates_2', [], 'rawspkrates_2_drug', [],... %<- ???
    'rate_resmpl_2', [], 'rate_resmpl_2_drug', [],... <- resampled firing rates
    'ratepar_2', [], 'ratepar_2_drug', [], ... %<- corresponding stimulus parameters
    'nrep_2', [], 'nrep_2_drug', [], ... %<- number of repeats per stimulus condition (only rewarded trials)
    'tcdiff', [], 'tcdiff_drug', [], ... %<- difference in tuning curve
    'tcdiff_1', [], 'tcdiff_1_drug', [], ... %<- difference in tuning curve
    'tcdiff_2', [], 'tcdiff_2_drug', [], ... %<- difference in tuning curve
    'ff', [], 'ff_drug', [], ... %<- struct for fano factor analysis
    'ff_1', [], 'ff_1_drug', [], ... %<- struct for fano factor analysis
    'ff_2', [], 'ff_2_drug', [], ... %<- struct for fano factor analysis
    'ppsz_mu', [],'ppsz_mu_drug', [], ... %<- fiel for pupil size analysis
    'resvars', [], 'resvars_drug', [], ... %<- ???
    'sdfs', {}, 'sdfs_drug', {}, ...%<- ???
    'times', [], 'times_drug', [], ...%<- ???
    'resvars_1', [], 'resvars_1_drug', [], ... %<- ???
    'sdfs_1', {}, 'sdfs_1_drug', {}, ...%<- ???
    'times_1', [], 'times_1_drug', [], ...%<- ???
    'resvars_2', [], 'resvars_2_drug', [], ... %<- ???
    'sdfs_2', {}, 'sdfs_2_drug', {}, ...%<- ???
    'times_2', [], 'times_2_drug', [], ...%<- ???
    'ismango', [], ...%<- boolean value indicating the monkey
    'others', [], ... %<- ???
    'ed', [], ... %<- electrode depth
    'bridx', [], ... %<- ???
    'isi_frct', [],... %<- ???
    'x0', [], 'y0', [], ... %<- stimulus position
    'ecc', [],... %<- ???
    'pfi', [], 'pfi_drug', [], ... %<- index for the preferred stimulus (highest response)
    'upfi', [], 'upfi_drug', [], ...%<- index for the unpreferred stimulus (lowest response)
    'tf_f1f0', [], ... %<- ???
    'p_anova', [], 'p_anova_drug', [], ... %<- ANOVA for significant tuning
    'p_anova_1', [], 'p_anova_1_drug', [], ... %<- ANOVA for significant tuning
    'p_anova_2', [], 'p_anova_2_drug', [], ... %<- ANOVA for significant tuning
    'phasesel', [], 'phasesel_drug', [], ... %<- ???
    'phasesel_1', [], 'phasesel_1_drug', [], ... %<- ???
    'phasesel_2', [], 'phasesel_2_drug', [], ... %<- ???
    'expduration', [], 'expduration_drug', [], ... %<- ???
    'reg_bootstrp', [],... %<- ???
    'psth', {}, 'psth_drug', {}, ... %<- matrix containing the smoothed PSTH
    'electrodebroken', [],... %<- boolean value. true if the electrode was broken
    'electrodebroken_excl', [], ... %<- boolean value. true if the electrode was broken and should be excluded from the analysis
    'electrodebroken_incl_underrest', [], ...%<- boolean value. true if the electrode was broken but should be included to the analysis
    'mntc_CI_base', [], 'mntc_CI_drug', [],...
    'c0rate', [], 'c0rate_drug', [], ... %<- ???
    'c0geomn', [], 'c0geomn_drug', [], ...%<- ???
    'c0geomn_2nd', [], 'c0geomn_2nd_drug', [], ... %<- ???
    'trials_c0', [], 'trials_c0_drug', [],... %<- ???
    'trials_c1', [], 'trials_c1_drug', [], ....%<- ???
    'is5HT', [], ... %<- boolean value, true when 5HT was applied (otherwise NaCl was applied)
    'isadapt', [], ... %<- true if the experiment was with an adaptation period
    'isRC', [], ... %<- is a flashed grating experiment, demands reverse corr. subspace analysis
    'isc2', [], ... %<- true if it's a cluster 2 spike sorted file
    'cluster', [], ... %<- more specifically
    'rsc', [], 'rsc_drug', [], ... %<- Pearson's correlation coefficient for z-scored spike counts of single and multiunit activity
    'prsc', [], 'prsc_drug', [], ... %<- corresponding p values
    'rsc_2nd', [], 'rsc_2nd_drug', [], ... %<- same computation for the second half odf the experiments onlz
    'prsc_2nd', [], 'prsc_2nd_drug', [], ...
    'rsig', [], 'rsig_drug', [], ... %<- Pearson's correlation coefficient for raw tuning curves of single and multiunit activity (signal correlation)
    'prsig', [], 'prsig_drug', [], ... %<- corresponding p values
    'rsig_2nd', [], 'rsig_2nd_drug', [], ... %<- same computation for the second half odf the experiments onlz
    'prsig_2nd', [], 'prsig_2nd_drug', [], ... 
    'rsc_1', [], 'rsc_1_drug', [], ... %<- Pearson's correlation coefficient for z-scored spike counts of single and multiunit activity
    'prsc_1', [], 'prsc_1_drug', [], ... %<- corresponding p values
    'rsc_2nd_1', [], 'rsc_2nd_1_drug', [], ... %<- same computation for the second half odf the experiments onlz
    'prsc_2nd_1', [], 'prsc_2nd_1_drug', [], ...
    'rsig_1', [], 'rsig_1_drug', [], ... %<- Pearson's correlation coefficient for raw tuning curves of single and multiunit activity (signal correlation)
    'prsig_1', [], 'prsig_1_drug', [], ... %<- corresponding p values
    'rsig_2nd_1', [], 'rsig_2nd_1_drug', [], ... %<- same computation for the second half odf the experiments onlz
    'prsig_2nd_1', [], 'prsig_2nd_1_drug', [], ... 
    'rsc_2', [], 'rsc_2_drug', [], ... %<- Pearson's correlation coefficient for z-scored spike counts of single and multiunit activity
    'prsc_2', [], 'prsc_2_drug', [], ... %<- corresponding p values
    'rsc_2nd_2', [], 'rsc_2nd_2_drug', [], ... %<- same computation for the second half odf the experiments onlz
    'prsc_2nd_2', [], 'prsc_2nd_2_drug', [], ...
    'rsig_2', [], 'rsig_2_drug', [], ... %<- Pearson's correlation coefficient for raw tuning curves of single and multiunit activity (signal correlation)
    'prsig_2', [], 'prsig_2_drug', [], ... %<- corresponding p values
    'rsig_2nd_2', [], 'rsig_2nd_2_drug', [], ... %<- same computation for the second half odf the experiments onlz
    'prsig_2nd_2', [], 'prsig_2nd_2_drug', [], ... 
    'wdt', [], ... %<- width of the spike wave form
    'figname',  [], ... %<- generic figure name, indicating original experiments
    'figpath', [], ... %<- ???
    'fig_tc', [], ... %<- filename for the tuning curve plot
    'fig_regl', [], ... %<-  filename for the type-II regression plot
    'fig_raster', [], ... %<- filename for the raster plot
    'fig_waveform', [], ... %<- filename for the spike waveform
    'fig_sdfs', [], ... %<- filename for figure showing the spike density function
    'fig_lfpPow', [], ... %<- dummy filename for LFP power spectrum
    'fig_lfpSpec', [], ... %<- dummy filename for LFP spectrogram
    'fig_lfpAvgSpk', [], ... %<- dummy filename for spike triggered average LFP
    'fig_bri', [], ... % filename for the bursting analysis
    'fig_ppszTraj', [], ... %<- filename for the figure showing the pupil size trajectory
    'fig_ppszHist', [], ... %<- filename for the figure showing the pupil size histogram
    'fig_phase', [], ... %<- figure filename for the results of the f1/f0 analysis for the experiments without tf manipulation
    'fig_phasetf', [], ... %<- figure filename for the results of the TF experiment
    'fig_varxtime', [], ... %<- right now, only a dummy 
    'fig_noisecorr', [], ... %<- figure filename for noise correlation plots
    'fig_recovery', [], ...%<- figure filename that shows the average activity of a unit for each experimetn within the session
    'fig_intr',[], ... % <- figure filename that shows interaction between pupil and drug on tuning properties
    'fig_eye', [], ... %<- figure filename that shows fixation precision and microsaccades
    'stimdur', [], ... %<- the stimulus presentation duration of the RC flashed grating (10 or 30ms)
    'tf', [], ... %<- the temporal frequency of the stimulus
    'valid', [], ... %<- experiments with the best type-II regression fit within a stimulus dimension of a unit
    'isdominant', [], ... %<- is the experiment with stimuli shown to the dominant eye
    'cmpExp', [], ... %<- ???
    'RFwx', [], ... %<- estimated receptive field width in x dimension
    'RFwy', [], ... %<- estimated receptive field width in y dimension
    'RFw', [], ... %<- average receptive field width in 2D
    'RFw_corr', [], ... %<- eccentricity corrected receptive field width
    'spkqual_base', [], 'spkqual_drug', [], ... %<- spike sorting quality, derived from the google spreadsheet
    'gaussr2', [], 'gaussr2_drug', [], ... %<- tuning fit quality as measured by the explained variance
    'dn', [], ... %<- preceeding number of drug experiment within the unit recording
    'dt_cum', [], ... %<- cummulated time of 5HT application before this experiment for one unit
    'dn_id', [], ... 
    'dn_session', [],... %<- preceeding number of drug experiment in a recording session
    'dt_id', [],...
    'dt_cum_id', [],...
    'dt_cum_session', [],...
    'gridX', [], 'gridY', [], ... %<- grid position - note that the information in the ex file are not always correct. look uo the google spreadsheet for the correct position.
    'expstrt', [], ... %<- time stamp of the first trial start
    'ret2base', [], ... %<- boolean value, true if this baseline experiment recovered from the effect of preceeding drug application 
    'iscmp2base', [], ... %<- ???
    'reg_slope', [], ... %<- result of the latency control analysis using subsampling. Slope of the regression fit Latency ~ Noise (Baseline SD in SDFs)
    'reg_off', [], ... %<- corresponding offset of the regression
    'noise', [], 'noise_drug', [], ... %<- noise (SD) across SDFs, RC experiments only
    'lat_c', [], 'lat_drug_c', [], ... %<- noise corrected latency estimate, RC experiments only
    'predicted_lat_change', [], ... %<- predicted latency change based on the change in baseline noise
    'fig_latjackknife', [], ... %<- figure of the experiment subsampling for the latency control analysis
    'fixationspan', [], 'fixationspan_drug', [], ... %<- fixation precision, theta, symmetry index
    'variance_eye', [], 'variance_eye_drug', [],...%<- fixation precision as variance x and y
    'recentering_win', [], 'recentering_win_drug', [], ... % <- window position due to recentering
    'recentering_eye', [], 'recentering_eye_drug', [], ... % <- eye position due to recentering
    'microsac_amplitude', [], 'microsac_amplitude_drug', [], ...% <- microsaccade amplitude
    'microsac_peakv', [], 'microsac_peakv_drug', [], ... %<- microsaccade peak velocity
    'microsac_duration', [], 'microsac_duration_drug', [], ... % <- microsaccade duration
    'microsac_angle', [], 'microsac_angle_drug', [], ... % <- microsaccade angle
    'microsac_counts', [], 'microsac_counts_drug', [], ... % <- microsaccade counts
    'fixationspan_1', [], 'fixationspan_1_drug', [], ... %<- fixation precision, theta, symmetry index
    'variance_eye_1', [], 'variance_eye_1_drug', [],...%<- fixation precision as variance x and y
    'recentering_win_1', [], 'recentering_win_1_drug', [], ... % <- window position due to recentering
    'recentering_eye_1', [], 'recentering_eye_1_drug', [], ... % <- eye position due to recentering
    'microsac_amplitude_1', [], 'microsac_amplitude_1_drug', [], ...% <- microsaccade amplitude
    'microsac_peakv_1', [], 'microsac_peakv_1_drug', [], ... %<- microsaccade peak velocity
    'microsac_duration_1', [], 'microsac_duration_1_drug', [], ... % <- microsaccade duration
    'microsac_angle_1', [], 'microsac_angle_1_drug', [], ... % <- microsaccade angle
    'microsac_counts_1', [], 'microsac_counts_1_drug', [], ... % <- microsaccade counts
    'fixationspan_2', [], 'fixationspan_2_drug', [], ... %<- fixation precision, theta, symmetry index
     'variance_eye_2', [], 'variance_eye_2_drug', [],...%<- fixation precision as variance x and y
    'recentering_win_2', [], 'recentering_win_2_drug', [], ... % <- window position due to recentering
    'recentering_eye_2', [], 'recentering_eye_2_drug', [], ... % <- eye position due to recentering
    'microsac_amplitude_2', [], 'microsac_amplitude_2_drug', [], ...% <- microsaccade amplitude
    'microsac_peakv_2', [], 'microsac_peakv_2_drug', [], ... %<- microsaccade peak velocity
    'microsac_duration_2', [], 'microsac_duration_2_drug', [], ... % <- microsaccade duration
    'microsac_angle_2', [], 'microsac_angle_2_drug', [], ... % <- microsaccade angle
    'microsac_counts_2', [], 'microsac_counts_2_drug', [] ... % <- microsaccade counts
    ); 
   


%% parse input arguments
j = 1;
while j<=length(varargin)
    switch varargin{j}
        case 'fname'
            fname = varargin{j+1};
        case 'exinfo'
            exinfo = varargin{j+1};
        case 'figdir'
            figdir.pre = varargin{j+1};
    end
    j = j+2;
end


%% predefine the figure folder directories

figdir.TC = [figdir.pre 'TC'];
figdir.Regression = [figdir.pre 'TypeIIRegression'];
figdir.Raster = [figdir.pre 'Rasterplot'];
figdir.Waveform = [figdir.pre 'SpikeWaveform'];
figdir.PPSZdev = [figdir.pre 'PupilSizeDev'];
figdir.PPSZhist = [figdir.pre 'PupilSizeHist'];
figdir.RC = [figdir.pre 'RevCor_plot'];
figdir.lfp = [figdir.pre 'LFP'];
figdir.BRI = [figdir.pre 'BRI'];
figdir.Phase = [figdir.pre 'PhaseSelectivity'];
figdir.PSTH = [figdir.pre 'PSTH'];
figdir.Variability = [figdir.pre 'Variability'];
figdir.Interaction = [figdir.pre 'Interaction'];
figdir.Eye = [figdir.pre 'Eye'];


checkdir(figdir); % check whether the folders already exist and, if not, create the folders


%% handle broken electrodes
% all unit IDs recorded with a broken electrode 
broken_id = [12.5  13.5  25.5  28.5  29.5  30.5  44.5  45.5  46.5  47.5  ...
    48.5  49.5  51.5  59.5  60.5  66.5  67.5  68.5  69.5  70.5  71.5  ...
    78.5  80.5  81.5  82.5  91.5  92.5  93.5  98.5  99.5  100.5  101.5 ...
    102.5  103.5  111.5  112.5  113.5  116.5  117.5  118.5  119.5  125.5 ...
    126.5  135.5  136.5  137.5  138.5  139.5  144.5  145.5  146.5  147.5  ...
    148.5  151.5  152.5  155.5  160.5  161.5  162.5  163.5  164.5  165.5  ...
    166.5  167.5  168.5  169.5  174.5  175.5  176.5  177.5  189.5  190.5  ...
    191.5  192.5  193.5  194.5  195.5  196.5  197.5  198.5  199.5  200.5  ...
    205.5 207.5 211.5 213.5 219.5 237.5 238.5 239.5:241.5 253.5 254.5 257.5 ...
    64  65  66  115  116  124  138  139  142  143  144  145  158  159 ...
    160 161  168  169  179  180  181 182  183  186  187  188  196 ...
    197 198  199  206 208  214  235  244  254  255 259  260  265 ...
    266 267  268  269  272  275  281  291  292  293  308 309  312  313 ...
    314  315  326  327 ];
%all units recorded with a broken electrode on the 5HT side which are to
%be excluded from the analysis
broken_excl_5HT = [115  116  187  196  197  198  214 281  308  309 ...
    [44, 71, 78, 11:13 66:70, 80:82, 101:103, 107:112, 135, 151:152, 160 161, ...
    175:177, 189:200, 212:216,  239:241 253 254 257]+0.5];
%all units recorded with a broken electrode on the NaCl side which are to
%be excluded from the analysis
broken_excl_nacl = [115  116  187  198  199  206  208 214  235 ...
    244  259  260  265  266  267  268  269 272  275  281  308  309 ...
    [11:13, 24, 30, 45:49, 59, 60, 66:70,71, 80:82, 91:93, 98:103, 107:112, ...
    113, 116:119, 135, 136:139, 144:148, 151, 152, 160:164, 175:177, 189:200, ...
    212:216, 237 238 239:241]+0.5];
%all units recorded with a broken electrode on the 5HT side which are to
%be included to the analysis
broken_incl_5HT_rest = broken_id(~ismember(broken_excl_5HT, broken_id));


%% load file directories
% open file
fileID = fopen(fname, 'r');
SU_dir = textscan(fileID, '%s'); F = SU_dir{1};
fclose(fileID);

% divide the filenames into baseline and drug experiments
idx = find( ~contains(F, '5HT') & ~contains(F, 'NaCl') );  % indices of all baseline experiments
F_base = F(idx); % filenames of the baseline experiments
F_drug = F(idx+1); % filenames of the drug (5HT or NaCl) experiments


%% assign experiment information to the exinfo struct


% loop through each experiment pair and derive information about the
% identifiers, e.g. ID, dose, etc. but also assign the figure directories
% and the results of the electrode inclusion criteria
for i = 1:length(F_drug)
    
    fprintf('WORKING ON ROW %1i \n', i);
    %----------------------------------------------------------- load files
    
    fname = F_base{i};
    fname_drug = F_drug{i};
    
    % load c1/2 and c0 ex files
    ex0 = loadCluster(fname, 'reward', false, 'loadlfp', false);    
    ex2 = loadCluster(fname_drug,'reward', false, 'loadlfp', false);
    
    
    %------------------------------- derive information from the filename
    % ...the applied drug
    if contains(fname_drug, 'NaCl')
        drugname = 'NaCl';
    elseif contains(fname_drug, '5HT')
        drugname = '5HT';
    else
       error('Could not determine the applied drug.');
    end
    
    % ...the cluster
    if contains(fname_drug, 'c2')
        spkcluster = 'c2';
    elseif contains(fname_drug, 'c1')
        spkcluster = 'c1';
    else
       error('Could not determine the spike cluster.'); 
    end
    
    % ...the animal
    if contains(fname, 'mango')
        idxa = 1;
        monkey = 'ma';
        id_off = 0;
    elseif contains(fname, 'kaki')
        idxa = 0;
        monkey = 'ka';
        id_off = 0.5;
    else
       error('Could not determine the animal name.');
    end
    
    id = str2double(fname(15+idxa:17+idxa))+id_off;
    
    %--------------------------------- derive information from the ex files
    % ...electrode depth
    if isfield(ex0.Trials, 'ed')
        ed = ex0.Trials(1).ed;
    else
        ed = nan;
    end
    
    % stimulus types (usually one or two, e.g. or and me, or and co,...)
    stimtypes = {ex0.exp.e1.type, ex0.exp.e2.type};
    
    % stimulus position in x and y, and resulting eccentricity 
    x0 = ex0.stim.vals.x0;  
    y0 = ex0.stim.vals.y0;  
    ecc = sqrt( x0^2 + y0^2 ); 
    
    % applied drug dose, in nA
    [dose, dosernd] = getDose(ex2); 
    
    % recorded voltage
    volt = getVolt(ex2);    
    
    % electrode position !inconsistant with the google spreadsheet table
    [gridX, gridY, ] = getGridPosition( ex0 ); 
    
    % beginning of the experiment
    expStrt = min( [ex0.Trials(1).TrialStart, ex2.Trials(1).TrialStart] ); 

    
    % generic figure name
    figname_gen = [fname((19:42)+idxa) '_' fname_drug((37:42)+idxa) ...
        '_' stimtypes{1} 'x' stimtypes{2} ];
    
    
    
    %--------------------------------- derive information from the ex files

    electrodebroken  = ismember(id, broken_id);

    % check whether electrodes were broken and excluded
    if strcmp(drugname, '5HT')
        electrodebroken_excl = ismember(id, broken_excl_5HT);
        electrodebroken_incl_underrest = ismember(id, broken_incl_5HT_rest);
    else
        electrodebroken_excl = ismember(id, broken_excl_nacl);
        electrodebroken_incl_underrest = 0; % NaCl experiments with broken electrode are never included
    end
    
    %------------------------------------------------------------ ocularity
    % in case the ocularity was a second variation of the stimulus, each
    % condition is handled individually and is assigned to a row in the
    % result structure
    for ocul = intersect( unique([ex0.Trials.me]), unique([ex2.Trials.me]) )
        
        kk = kk +1;
        
        %----------------------------------------------------- general info
        exinfo(kk).id = str2double(fname(14+idxa:17+idxa))+id_off;
        exinfo(kk).idi = idi;   % index indicating the same file
        exinfo(kk).electrodebroken = any(broken_id == exinfo(kk).id);
        exinfo(kk).monkey = monkey;
        exinfo(kk).ismango = strcmp(monkey, 'ma');
        exinfo(kk).dose = dose;
        exinfo(kk).dosernd = dosernd;
        exinfo(kk).ed = ed;
        exinfo(kk).volt = volt;
        exinfo(kk).resistance = volt/dose; 
        
        exinfo(kk).x0 = x0;                   
        exinfo(kk).y0 = y0;
        exinfo(kk).ecc = ecc;
        exinfo(kk).gridX = gridX;
        exinfo(kk).gridY = gridY;
        exinfo(kk).expstrt = expStrt;
        exinfo(kk).expduration = ex0.Trials(end).TrialEnd - ex0.Trials(1).TrialStart;
        exinfo(kk).expduration_drug = ex2.Trials(end).TrialEnd - ex2.Trials(1).TrialStart;
        
        exinfo(kk).fname = fname;             exinfo(kk).fname_drug = fname_drug;
        exinfo(kk).date = getExDate(ex0);     exinfo(kk).date_drug = getExDate(ex2);
        
        
        exinfo(kk).ocul = ocul;
        exinfo(kk).drugname = drugname;
        
        
        exinfo(kk).param1 = stimtypes{1};      exinfo(kk).param2 = stimtypes{2};
        exinfo(kk).tf = ex0.stim.vals.tf;
        exinfo(kk).is5HT      = contains(drugname, '5HT');
        exinfo(kk).isadapt    = contains(fname_drug, 'adapt');
        exinfo(kk).isRC       = contains(fname_drug, 'RC');
        exinfo(kk).isc2       = strcmp(spkcluster, 'c2');
        exinfo(kk).cluster    = spkcluster;
        
        exinfo(kk).wdt = -1;
        exinfo(kk).lat = -10;
        exinfo(kk).lat_drug = -10;
        exinfo(kk).gaussr2 = -10;
        exinfo(kk).gaussr2_drug = -10;
        

        
        if exinfo(kk).isRC
            if isfield(ex0.stim.vals,'RCperiod')
                exinfo(kk).stimdur = ex0.stim.vals.RCperiod;
            else
                exinfo(kk).stimdur = 1;
            end
            add_RC = 'xRC';
        else
            add_RC = '';
        end
        
        
        %----------------------------------------------------- figure names
        
        switch ocul
            case -1
                figname = [figname_gen add_RC '_left_' drugname];
            case -11
                figname = [figname_gen add_RC '_left_' drugname];
            case 0
                figname = [figname_gen add_RC '_both_' drugname];
            case 1
                figname = [figname_gen add_RC '_right_' drugname];
        end
        
        
        exinfo(kk).figname        = figname;
        exinfo(kk).figpath        = fullfile(figdir.TC, '\', [figname '_summary.fig']);
        exinfo(kk).fig_tc         = fullfile(figdir.TC, '\', [figname '_tc.fig']);
        exinfo(kk).fig_regl       = fullfile(figdir.Regression, '\', [figname '_regl.fig']);
        exinfo(kk).fig_raster     = fullfile(figdir.Raster, '\', [figname '_raster.fig']);
        exinfo(kk).fig_waveform   = fullfile(figdir.Waveform, '\', [figname '_waveform.fig']);
        
        exinfo(kk).fig_sdfs       = fullfile(figdir.RC, '\', [figname '_sdfs.fig']);
        
        exinfo(kk).fig_lfpPow     = fullfile(figdir.lfp, '\', [figname '_lfppow.fig']);
        exinfo(kk).fig_lfpSpec     = fullfile(figdir.lfp, '\', [figname '_lfpspec.fig']);
        exinfo(kk).fig_lfpAvgSpk   = fullfile(figdir.lfp, '\', [figname '_spktriglfp.fig']);
        
        
        exinfo(kk).fig_bri        = fullfile(figdir.BRI, '\', [figname '_bri.fig']);
        
        exinfo(kk).fig_ppszTraj   = fullfile(figdir.PPSZdev, '\', [figname '_devppsz_b4.fig']);
        exinfo(kk).fig_ppszHist   = fullfile(figdir.PPSZhist, '\', [figname '_histppsz.fig']);
        
        exinfo(kk).fig_phase   = fullfile(figdir.Phase, '\', [figname '_phase.fig']);
        exinfo(kk).fig_psth   = fullfile(figdir.PSTH, '\', [figname '_psth.fig']);
        
        exinfo(kk).fig_phasetf = fullfile(figdir.Phase, '\', [figname '_psth_tf.fig']);
        
        exinfo(kk).fig_varxtime = fullfile(figdir.Variability, '\', [figname '_varxtime.fig']);
        exinfo(kk).fig_noisecorr = fullfile(figdir.Variability, '\', [figname '_noisecorr.fig']);
        
        exinfo(kk).fig_intr = fullfile(figdir.Interaction, '\', [figname '_intr.fig']);
        exinfo(kk).fig_eye = fullfile(figdir.Eye, '\', [figname '_eye.fig']);
        
    end
    
    % assing whether the electrode were broken and included 
    exinfo(kk).electrodebroken  = electrodebroken;
    exinfo(kk).electrodebroken_excl = electrodebroken_excl;
    exinfo(kk).electrodebroken_incl_underrest = electrodebroken_incl_underrest;
    
    idi = idi+1;
    
end




% retrieve information from other sources

% phase selectivity --> TF experiment
exinfo = setPhaseSelTFexp( exinfo );

% RF size --> XPos / YPos experiment
exinfo = setReceptiveFieldSize( exinfo , 1 );

% experiment number in the recording session, spike sporting isolation quality, etc
% --> google spreadsheet
exinfo = addSortingValue(exinfo);
exinfo = addNumInExp(exinfo);

save(fullfile(figdir.Data, 'empty_exinfo.mat'), 'exinfo')

end




function checkdir(fdir)
% create the file directory, if it does not exist yet

fnames = fieldnames(fdir);
for i = 1:length(fnames)
    if ~exist(fdir.(fnames{i}), 'dir')
        mkdir(fdir.(fnames{i}));
    end
end
end




