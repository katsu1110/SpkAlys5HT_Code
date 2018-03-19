function exinfo = runExinfoAnalysis( varargin )
%exinfo = runExinfoAnalysis()
%
% Batch file calls all functions that perform the analytical computation on
% exinfo. exinfo is the result structure created by initExinfo(). If is is
% not given as input, initExinfo() is called. 
% 
%
% Optional arguments are:
% 'plot'            -  plots results and saves figures according to the
%                        destination in exinfo (specified in initExinfo)
% 'row_i_strt', XX  - start the for-loop at the row XX, helpful for
%                       debugging
% 'exinfo', XX      - where exinfo is the result structure. This avoids
%                       running initExinfo. 
% 'save'            - this input causes the output to be saved as
%                       exinfo.mat
%
% 
% 
% 
% This function calls the following functions forwarding any input argument
% given to runExinfoAnalysis
%   initExinfo      - initializes exinfo if none is given as input.
%   evalSingleEx    - returns single unit analysis results
%   evalBothEx      - returns the comparative analysis results
%
%     Arguments passed on to evalSingleDG:
%     'rsc'             - run the analysis of co-variability.
%     'ff'              - the analysis of variability is performed. 
%     'anova'           - run an ANOVA on the raw tuning curve to check the
%                           selectivity.
%     'tcfit'           - fit the descriptive functions to the tuning curve.
%     'tcheight'        - compute the amplitude of the tuning curve divided by
%                           the mean height.
%     'phasesel'        - run the computation of phase selectivity on this
%                           experiment.
% 
%     Arguments passed on to evalBothEx:
%       't2reg'         - fit the type-II linear regression
%       'wavewdt'       - compute the spike waveform's wave width
%      
% 
% 'all'         - performs all of these computations.
% 
%     Arguments passed on to evalSingleRC-RCsubspace:
%       'plot'                  - plot the sdfs and latency estimates
%       'lat_flag', boolean     - compute the ML latency estimate (default) or not.
%       'bootstrap_n', integer  - specify the number of bootstrap samples
%
% 
% The function further calls these functions  
%   getValidField           - determines the experiment with best
%                               type-II regression fit for experiments with
%                               the same stimulus type for the same unit
%                               cluster
%   getDominantEyeField     - determines the eye preference, important for
%                               experiments with varying ocularity
%   addStruct               - obsolete but important for the gui. adds the
%                               fields gaussr2(_drug) and sets the latency
%                               field to -10 if it is empty.
% 
% 
% @CL 


clc;
close all;
%% define default variables and parse input

% adding folder to the MATLAB search path takes some time
% I normally just do this ones and then comment it
disp('wait...')

addpath(genpath(pwd)); 
addpath(genpath('C:\Users\katsuhisa\Documents\code\analysis\integrated\interaction_project\GenAlyz_Code'));
addpath(genpath('C:\Users\katsuhisa\Documents\code\analysis\integrated\corinnas\HNetc'));

% addpath(genpath('Z:\Corinna\SharedCode\File Exchange Code')); % add all subfolders to the path

disp('all the paths added.')

rng(2384569);       % set seed for the random number generator. do not change this number. 
exinfo = [];

p_flag = false;     % not plotting the results is the default setting
i_strt = 1;         % starting index for the loop through exinfo
save_flag = false;   % not saving the result structure is default

j = 1;              
while  j<= length(varargin)
    switch varargin{j}
        case 'plot' % input string to call for plots
            p_flag = true;
            j = j+1;
            close all
        case 'row_i_strt' % input pair to start the loop at a specific row
            i_strt = varargin{j+1};
            j = j+2;
        case 'exinfo' % result structure. prevents initExinfo
            exinfo = varargin{j+1};
            varargin = varargin([1:j-1,j+2:end]);
        case 'save' % input string to save the result structure
            save_flag = true;
            j = j+1;
        otherwise
            j = j+1;
    end
end


% initialize the no exinfo structure in case there is none given as input
if isempty(exinfo)
    exinfo = initExinfo(varargin{:});
end

%% loop through each file and add other information
for kk = i_strt:length(exinfo)

%     if exinfo(kk).isRC % helpfull for debugging
%         continue;
%     end
    
    fprintf('WORKING ON ROW %1i, file %1.1f \n', kk, exinfo(kk).id);

    %--------------------------------------- operations on single ex files
    %baseline
    [ex0, args0] = evalSingleEx(exinfo(kk), exinfo(kk).fname, varargin{:});
    exinfo(kk) = assignInfo(exinfo(kk), '', args0{:});
    
    %drug
    [ex2, args2] = evalSingleEx(exinfo(kk), exinfo(kk).fname_drug, varargin{:});
    exinfo(kk) = assignInfo(exinfo(kk), '_drug', args2{:});
    
    
    %------------------------------------------ operations on both ex files
    exinfo(kk) = evalBothEx(ex0, ex2,  exinfo(kk), p_flag, varargin{:});
  
       
    %--------------------------------------------------------- plot results
    if any(strcmp(varargin, 'pupil')) || any(strcmp(varargin, 'all'))
        try
            % interaction btw. 5HT & pupil
            plotIntr(exinfo(kk));
        catch
            disp('error in plotIntr')
        end
    end
    
    if ~exinfo(kk).isRC && p_flag

        if ~exinfo(kk).isadapt
                
                % PSTH
                exinfo(kk) = psthPlot(exinfo(kk), ex0, ex2);
                
                % tuning curve
                tuningCurvePlot(exinfo(kk));        
                
%                 % interaction btw. 5HT & pupil
%                 plotIntr(exinfo(kk));
                
                % phase
                exinfo(kk) = phasePlot(exinfo(kk), ex0, ex2);            
            
           % these functions resort on variability and co-variability related fields of exinfo
           % to this end they are only performed if the analysis was
           % computed
           if any(strcmp(varargin, 'all')) ||...
                   (any(strcmp(varargin, 'ff')) && any(strcmp(varargin, 'rsc' )))
                
                % raster
                rasterPlot( exinfo(kk), ex0, ex2);
                
                % normalized FR
                znormplot(ex0, ex2, exinfo(kk));
                
                % FF
                VariabilityPlot(exinfo(kk), ex0, ex2);
           end
           
        end

        exinfo(kk) = getISI_All(exinfo(kk), ex0, ex2, p_flag);

    elseif exinfo(kk).isRC && p_flag
        rcPlot(exinfo(kk));
        tuningCurvePlot(exinfo(kk));  
    end
    
    %-------------------------------------- temp save, useful for debugging
%     if save_flag && mod(kk, 30)==0 
%         save('exinfo.mat', 'exinfo', '-v7.3'); 
%     end
end



%--------------------------------------------------------------- add fields
exinfo = getValidField(exinfo);
exinfo = getDominantEyeField(exinfo);

exinfo = addStruct(exinfo);

try
    [exinfo] = addAll(exinfo, 0);
    disp('new analysis on isolation quality and variability is on!')
catch
    disp('addAll had an error to be fixed...')
end

% save the result structure in a superordinate folder named Data
if save_flag
%     save('D:\Users\kk\interaction_project\dataset\Data\exinfo.mat', 'exinfo', '-v7.3'); 
    save('Z:\Katsuhisa\serotonin_project\dataset\Data\exinfo.mat', 'exinfo', '-v7.3'); 
end

% % store figures for fixation precision
% storeEye(exinfo)

end

function exinfo = assignInfo(exinfo, apx, varargin)
% arguments from evalSingleEx assigned to exinfo
%
% for this function it is important that input and output exinfo have the
% exact same fields.


j = 1;

while j<length(varargin)

    switch varargin{j}
        
        case 'lat'
            eval([ 'exinfo.lat' apx ' = varargin{j+1};']);
        case 'lat2Hmax'
            eval([ 'exinfo.lat2Hmax' apx ' = varargin{j+1};']);
        case 'fitparam'
            eval([ 'exinfo.fitparam' apx ' = varargin{j+1};']);
        case 'fitparam_1'
            eval([ 'exinfo.fitparam_1' apx ' = varargin{j+1};']);
        case 'fitparam_2'
            eval([ 'exinfo.fitparam_2' apx ' = varargin{j+1};']);
        case 'rateMN'
            eval([ 'exinfo.ratemn' apx ' = varargin{j+1};']);
        case 'rateMN_1'
            eval([ 'exinfo.ratemn_1' apx ' = varargin{j+1};']);
        case 'rateMN_2'
            eval([ 'exinfo.ratemn_2' apx ' = varargin{j+1};']);
        case 'rateVARS'
            eval([ 'exinfo.ratevars' apx ' = varargin{j+1};']);
        case 'rateVARS_1'
            eval([ 'exinfo.ratevars_1' apx ' = varargin{j+1};']);
        case 'rateVARS_2'
            eval([ 'exinfo.ratevars_2' apx ' = varargin{j+1};']);
        case 'ratePAR'
            eval([ 'exinfo.ratepar' apx ' = varargin{j+1};']);
        case 'ratePAR_1'
            eval([ 'exinfo.ratepar_1' apx ' = varargin{j+1};']);
        case 'ratePAR_2'
            eval([ 'exinfo.ratepar_2' apx ' = varargin{j+1};']);
        case 'rateSME'
            eval([ 'exinfo.ratesme' apx ' = varargin{j+1};']);
        case 'rateSME_1'
            eval([ 'exinfo.ratesme_1' apx ' = varargin{j+1};']);
        case 'rateSME_2'
            eval([ 'exinfo.ratesme_2' apx ' = varargin{j+1};']);
        case 'rateSD'
            eval([ 'exinfo.ratesd' apx ' = varargin{j+1};']);
        case 'rateSD_1'
            eval([ 'exinfo.ratesd_1' apx ' = varargin{j+1};']);
        case 'rateSD_2'
            eval([ 'exinfo.ratesd_2' apx ' = varargin{j+1};']);
        case 'rawspkrates'
            eval([ 'exinfo.rawspkrates' apx ' = varargin{j+1};']);
        case 'rawspkrates_1'
            eval([ 'exinfo.rawspkrates_1' apx ' = varargin{j+1};']);
        case 'rawspkrates_2'
            eval([ 'exinfo.rawspkrates_2' apx ' = varargin{j+1};']);
        case 'rate_resmpl'
            eval([ 'exinfo.rate_resmpl' apx ' = varargin{j+1};']);
        case 'rate_resmpl_1'
            eval([ 'exinfo.rate_resmpl_1' apx ' = varargin{j+1};']);
        case 'rate_resmpl_2'
            eval([ 'exinfo.rate_resmpl_2' apx ' = varargin{j+1};']);
        case 'ff'
            eval([ 'exinfo.ff' apx ' = varargin{j+1};']);
        case 'ff_1'
            eval([ 'exinfo.ff_1' apx ' = varargin{j+1};']);
        case 'ff_2'
            eval([ 'exinfo.ff_2' apx ' = varargin{j+1};']);
        case 'tcdiff'
            eval([ 'exinfo.tcdiff' apx ' = varargin{j+1};']);
        case 'tcdiff_1'
            eval([ 'exinfo.tcdiff_1' apx ' = varargin{j+1};']);
        case 'tcdiff_2'
            eval([ 'exinfo.tcdiff_2' apx ' = varargin{j+1};']);
        case 'rsc'
            eval([ 'exinfo.rsc' apx ' = varargin{j+1};']);
        case 'rsc_1'
            eval([ 'exinfo.rsc_1' apx ' = varargin{j+1};']);
        case 'rsc_2'
            eval([ 'exinfo.rsc_2' apx ' = varargin{j+1};']);
        case 'prsc'
            eval([ 'exinfo.prsc' apx ' = varargin{j+1};']);
        case 'prsc_1'
            eval([ 'exinfo.prsc_1' apx ' = varargin{j+1};']);
        case 'prsc_2'
            eval([ 'exinfo.prsc_2' apx ' = varargin{j+1};']);
        case 'rsig'
            eval([ 'exinfo.rsig' apx ' = varargin{j+1};']);
        case 'rsig_1'
            eval([ 'exinfo.rsig_1' apx ' = varargin{j+1};']);
        case 'rsig_2'
            eval([ 'exinfo.rsig_2' apx ' = varargin{j+1};']);
        case 'prsig'
            eval([ 'exinfo.prsig' apx ' = varargin{j+1};']);
        case 'prsig_1'
            eval([ 'exinfo.prsig_1' apx ' = varargin{j+1};']);
        case 'prsig_2'
            eval([ 'exinfo.prsig_2' apx ' = varargin{j+1};']);
        case 'rsc_2nd'
            eval([ 'exinfo.rsc_2nd' apx ' = varargin{j+1};']);
        case 'rsc_2nd_1'
            eval([ 'exinfo.rsc_2nd_1' apx ' = varargin{j+1};']);
        case 'rsc_2nd_2'
            eval([ 'exinfo.rsc_2nd_2' apx ' = varargin{j+1};']);
        case 'prsc_2nd'
            eval([ 'exinfo.prsc_2nd' apx ' = varargin{j+1};']);
        case 'prsc_2nd_1'
            eval([ 'exinfo.prsc_2nd_1' apx ' = varargin{j+1};']);
        case 'prsc_2nd_2'
            eval([ 'exinfo.prsc_2nd_2' apx ' = varargin{j+1};']);
        case 'rsig_2nd'
            eval([ 'exinfo.rsig_2nd' apx ' = varargin{j+1};']);
        case 'rsig_2nd_1'
            eval([ 'exinfo.rsig_2nd_1' apx ' = varargin{j+1};']);
        case 'rsig_2nd_2'
            eval([ 'exinfo.rsig_2nd_2' apx ' = varargin{j+1};']);
        case 'prsig_2nd'
            eval([ 'exinfo.prsig_2nd' apx ' = varargin{j+1};']);
        case 'prsig_2nd_1'
            eval([ 'exinfo.prsig_2nd_1' apx ' = varargin{j+1};']);
        case 'prsig_2nd_2'
            eval([ 'exinfo.prsig_2nd_2' apx ' = varargin{j+1};']);
        case 'resvars'
            eval([ 'exinfo.resvars' apx ' = varargin{j+1};']);
        case 'expduration'
            eval([ 'exinfo.expduration' apx ' = varargin{j+1};']);
        case 'sdfs'
            eval([ 'exinfo.sdfs' apx ' = varargin{j+1};']);
        case 'times'
            eval([ 'exinfo.times' apx ' = varargin{j+1};']);
        case 'resdur'
            eval([ 'exinfo.dur' apx ' = varargin{j+1};']);
        case 'ed'
            exinfo.ed = varargin{j+1};
        case 'eX'
            exinfo.eX = varargin{j+1};
        case 'eY'
            exinfo.eY = varargin{j+1};
        case  'phasesel'
            eval([ 'exinfo.phasesel' apx ' = varargin{j+1};']);
        case  'phasesel_1'
            eval([ 'exinfo.phasesel_1' apx ' = varargin{j+1};']);
        case  'phasesel_2'
            eval([ 'exinfo.phasesel_2' apx ' = varargin{j+1};']);
        case  'psth'
            eval([ 'exinfo.psth' apx ' = varargin{j+1};']);
        case 'p_anova'
             eval([ 'exinfo.p_anova' apx ' = varargin{j+1};']);
        case 'p_anova_1'
             eval([ 'exinfo.p_anova_1' apx ' = varargin{j+1};']);
        case 'p_anova_2'
             eval([ 'exinfo.p_anova_2' apx ' = varargin{j+1};']);
        case 'nrep'
            eval([ 'exinfo.nrep' apx ' = varargin{j+1};']);
        case 'c0rate'
            eval([ 'exinfo.c0rate' apx ' = varargin{j+1};']);
        case 'c0geomn'
            eval([ 'exinfo.c0geomn' apx ' = varargin{j+1};']);
        case 'c0geomn_2nd'
            eval([ 'exinfo.c0geomn_2nd' apx ' = varargin{j+1};']);
        case 'trials_c0'
            eval([ 'exinfo.trials_c0' apx ' = varargin{j+1};']);
        case 'trials_c1'
            eval([ 'exinfo.trials_c1' apx ' = varargin{j+1};']);
        case 'fixationspan'
            eval([ 'exinfo.fixationspan' apx ' = varargin{j+1};']);
        case 'variance_eye'
            eval([ 'exinfo.variance_eye' apx ' = varargin{j+1};']);
        case 'recentering_win'
            eval(['exinfo.recentering_win' apx '= varargin{j+1};']);
       case 'recentering_eye'
            eval([ 'exinfo.recentering_eye' apx ' = varargin{j+1};']);
        case 'microsac_amplitude'
            eval([ 'exinfo.microsac_amplitude' apx ' = varargin{j+1};']);
        case 'microsac_peakv'
            eval([ 'exinfo.microsac_peakv' apx ' = varargin{j+1};']);
        case 'microsac_duration'
            eval([ 'exinfo.microsac_duration' apx ' = varargin{j+1};']);
        case 'microsac_angle'
            eval([ 'exinfo.microsac_angle' apx ' = varargin{j+1};']);
        case 'microsac_counts'
            eval(['exinfo.microsac_counts' apx '= varargin{j+1};']);
        case 'fixationspan_1'
            eval([ 'exinfo.fixationspan_1' apx ' = varargin{j+1};']);
        case 'variance_eye_1'
            eval([ 'exinfo.variance_eye_1' apx ' = varargin{j+1};']);
        case 'recentering_win_1'
            eval(['exinfo.recentering_win_1' apx '= varargin{j+1};']);
       case 'recentering_eye_1'
            eval([ 'exinfo.recentering_eye_1' apx ' = varargin{j+1};']);
        case 'microsac_amplitude_1'
            eval([ 'exinfo.microsac_amplitude_1' apx ' = varargin{j+1};']);
        case 'microsac_peakv_1'
            eval([ 'exinfo.microsac_peakv_1' apx ' = varargin{j+1};']);
        case 'microsac_duration_1'
            eval([ 'exinfo.microsac_duration_1' apx ' = varargin{j+1};']);
        case 'microsac_angle_1'
            eval([ 'exinfo.microsac_angle_1' apx ' = varargin{j+1};']);
        case 'microsac_counts_1'
            eval(['exinfo.microsac_counts_1' apx '= varargin{j+1};']);
        case 'fixationspan_2'
            eval([ 'exinfo.fixationspan_2' apx ' = varargin{j+1};']);
        case 'variance_eye_2'
            eval([ 'exinfo.variance_eye_2' apx ' = varargin{j+1};']);
        case 'recentering_win_2'
            eval(['exinfo.recentering_win_2' apx '= varargin{j+1};']);
       case 'recentering_eye_2'
            eval([ 'exinfo.recentering_eye_2' apx ' = varargin{j+1};']);
        case 'microsac_amplitude_2'
            eval([ 'exinfo.microsac_amplitude_2' apx ' = varargin{j+1};']);
        case 'microsac_peakv_2'
            eval([ 'exinfo.microsac_peakv_2' apx ' = varargin{j+1};']);
        case 'microsac_duration_2'
            eval([ 'exinfo.microsac_duration_2' apx ' = varargin{j+1};']);
        case 'microsac_angle_2'
            eval([ 'exinfo.microsac_angle_2' apx ' = varargin{j+1};']);
        case 'microsac_counts_2'
            eval(['exinfo.microsac_counts_2' apx '= varargin{j+1};']);

    end
    j = j+2;
end

end

