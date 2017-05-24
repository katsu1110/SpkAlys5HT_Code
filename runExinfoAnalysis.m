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
% This function calls the following functions forwarding any input argument
% given to runExinfoAnalysis
%   evalSingleEx    - returns single unit analysis results
%   evalBothEx      - returns the comparative analysis results
% 
% It further calls these functions  
%   getValidField           - ???
%   getDominantEyeField     - determines the eye preference, important for
%                               experiments with varying ocularity
% - setReceptiveFieldSize   - determines the RF width by finding and
%                               analyzing the YPos.mat and XPos.mat files
% - addSortingValue         - retrieves the spike sorting quality
% - addNumInExp             - retrieves information about the preceeding
%                               drug experiments
% - addStruct               - obsolete but important for the gui. adds the
%                               fields gaussr2(_drug) and sets the latency
%                               field to -10 if it is empty.
% 
% 
% @CL 



%% define default variables and parse input

exinfo = [];
rng(2384569);       % set seed for the random number generator. do not change this number. 
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
            j = j+2;
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
    exinfo(kk).tf = ex0.stim.vals.tf;
    exinfo(kk) = evalBothEx(ex0, ex2,  exinfo(kk), p_flag, varargin{:});
  
       
    %--------------------------------------------------------- plot results

    if ~exinfo(kk).isRC && p_flag

        if ~exinfo(kk).isadapt

           exinfo(kk) = new_psthPlot_red(exinfo(kk), ex0, ex2);
           rasterPlot( exinfo(kk), ex0, ex2);
           tuningCurvePlot(exinfo(kk));        
           znormplot(ex0, ex2, exinfo(kk));
           exinfo(kk) = phasePlot(exinfo(kk), ex0, ex2);
           
           VariabilityPlot(exinfo(kk), ex0, ex2);
        end

        exinfo(kk) = getISI_All(exinfo(kk), ex0, ex2, p_flag);

    elseif exinfo(kk).isRC && p_flag
        rcPlot(exinfo(kk));
        tuningCurvePlot(exinfo(kk));  
    end
    
    %-------------------------------------- temp save, useful for debugging
%     if save_flag && mod(kk, 30)==0 
%         save(['exinfo' fig_suffix '.mat'], 'exinfo', '-v7.3'); 
%     end
end

%--------------------------------------------------------------- add fields
exinfo = getValidField(exinfo);
exinfo = getDominantEyeField(exinfo);
exinfo = setReceptiveFieldSize( exinfo );
exinfo = addSortingValue(exinfo);
exinfo = addNumInExp(exinfo);
exinfo = addStruct(exinfo);


% save the result structure in a superordinate folder named Data
if save_flag
    cfolder = cd('..');
    datadir = fullfile(cd, 'Data\');%folder destination 
    cd(cfolder);
    save(fillfile(datadir, ['exinfo' fig_suffix '.mat']), 'exinfo', '-v7.3'); 
end

end

function exinfo = assignInfo(exinfo, apx, varargin)
% arguments from evalSingleEx assigned to exinfo
j = 1;

while j<length(varargin)

    switch varargin{j}
        
        case 'lat'
            eval([ 'exinfo.lat' apx ' = varargin{j+1};']);
        case 'lat2Hmax'
            eval([ 'exinfo.lat2Hmax' apx ' = varargin{j+1};']);
        case 'fitparam'
            eval([ 'exinfo.fitparam' apx ' = varargin{j+1};']);
        case 'rateMN'
            eval([ 'exinfo.ratemn' apx ' = varargin{j+1};']);
        case 'rateVARS'
            eval([ 'exinfo.ratevars' apx ' = varargin{j+1};']);
        case 'ratePAR'
            eval([ 'exinfo.ratepar' apx ' = varargin{j+1};']);
        case 'rateSME'
            eval([ 'exinfo.ratesme' apx ' = varargin{j+1};']);
        case 'rateSD'
            eval([ 'exinfo.ratesd' apx ' = varargin{j+1};']);
        case 'rawspkrates'
            eval([ 'exinfo.rawspkrates' apx ' = varargin{j+1};']);
        case 'rate_resmpl'
            eval([ 'exinfo.rate_resmpl' apx ' = varargin{j+1};']);
        case 'ff'
            eval([ 'exinfo.ff' apx ' = varargin{j+1};']);
        case 'tcdiff'
            eval([ 'exinfo.tcdiff' apx ' = varargin{j+1};']);
        case 'rsc'
            eval([ 'exinfo.rsc' apx ' = varargin{j+1};']);
        case 'prsc'
            eval([ 'exinfo.prsc' apx ' = varargin{j+1};']);
        case 'rsig'
            eval([ 'exinfo.rsig' apx ' = varargin{j+1};']);
        case 'prsig'
            eval([ 'exinfo.prsig' apx ' = varargin{j+1};']);
        case 'rsc_2nd'
            eval([ 'exinfo.rsc_2nd' apx ' = varargin{j+1};']);
        case 'prsc_2nd'
            eval([ 'exinfo.prsc_2nd' apx ' = varargin{j+1};']);
        case 'rsig_2nd'
            eval([ 'exinfo.rsig_2nd' apx ' = varargin{j+1};']);
        case 'prsig_2nd'
            eval([ 'exinfo.prsig_2nd' apx ' = varargin{j+1};']);
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
        case  'psth'
            eval([ 'exinfo.psth' apx ' = varargin{j+1};']);
        case 'p_anova'
             eval([ 'exinfo.p_anova' apx ' = varargin{j+1};']);
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
    end
    j = j+2;
end

end

