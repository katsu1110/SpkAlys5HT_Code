function exinfo = runExinfoAnalysis( varargin )
%exinfo = runExinfoAnalysis()
%
% 
%
%
%
% Optional arguments are:
% 'plot'            -  plots results and saves figures according to the
%                        destination in exinfo (specified in initExinfo)
% 'row_i_strt', XX  - start the for loop at the row XX, helpful for
%                       debugging
% 'exinfo', XX      - where exinfo is the result structure. This avoids
%                       running initExinfo. 
% 'save'            - this input causes the output to be saved as
%                       exinfo.mat
% 
% 
% 
% 
% @CL 



%% define default variables and parse input

exinfo = [];
rng(2384569);       % set seed for the random number generator. changing this number can cause deviating results of fitting algorithms. 
p_flag = false;     % no plotting is the default setting
i_strt = 1;         % row to start the loop through exinfo
save_flag = false;   % not saving the result is default

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


% in case there is no exinfo structure, initialize it 
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
    %base
    [ex0, args0] = evalSingleEx(exinfo(kk), exinfo(kk).fname, varargin{:});
    exinfo(kk) = assignInfo(exinfo(kk), '', args0{:});
    
    %drug
    [ex2, args2] = evalSingleEx(exinfo(kk), exinfo(kk).fname_drug, varargin{:});
    exinfo(kk) = assignInfo(exinfo(kk), '_drug', args2{:});

    %-------------------------------------- operations on both ex files
    exinfo(kk).tf = ex0.stim.vals.tf;
    exinfo(kk) = evalBothEx(ex0, ex2,  exinfo(kk), p_flag, varargin{:});
  
       
    %-------------------------------------- plot results

    if ~exinfo(kk).isRC

        if p_flag && ~exinfo(kk).isadapt

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
    
    %-------------------------------------- temp save
    if save_flag && mod(kk, 30)==0 
%         save(['exinfo' fig_suffix '.mat'], 'exinfo', '-v7.3'); 
    end
end

%------------------------------------------------------------- add fields
exinfo = getValidField(exinfo);
exinfo = getDominantEyeField(exinfo);
exinfo = setReceptiveFieldSize( exinfo );
exinfo = addSortingValue(exinfo);
exinfo = addNumInExp(exinfo);
exinfo = addStruct(exinfo);



if save_flag

    cfolder = cd('..');
    Datadir = fullfile(cd, 'Data\');%folder destination for data
    cd(cfolder);
    save(fillfile(Datadir, ['exinfo' fig_suffix '.mat']), 'exinfo', '-v7.3'); 

end

end

function info = assignInfo(info, apx, varargin)
% arguments assigned to in evalSingleEx
j = 1;

while j<length(varargin)

    switch varargin{j}
        
        case 'lat'
            eval([ 'info.lat' apx ' = varargin{j+1};']);
        case 'lat2Hmax'
            eval([ 'info.lat2Hmax' apx ' = varargin{j+1};']);
        case 'fitparam'
            eval([ 'info.fitparam' apx ' = varargin{j+1};']);
        case 'rateMN'
            eval([ 'info.ratemn' apx ' = varargin{j+1};']);
        case 'rateVARS'
            eval([ 'info.ratevars' apx ' = varargin{j+1};']);
        case 'ratePAR'
            eval([ 'info.ratepar' apx ' = varargin{j+1};']);
        case 'rateSME'
            eval([ 'info.ratesme' apx ' = varargin{j+1};']);
        case 'rateSD'
            eval([ 'info.ratesd' apx ' = varargin{j+1};']);
        case 'rawspkrates'
            eval([ 'info.rawspkrates' apx ' = varargin{j+1};']);
        case 'rate_resmpl'
            eval([ 'info.rate_resmpl' apx ' = varargin{j+1};']);
        case 'ff'
            eval([ 'info.ff' apx ' = varargin{j+1};']);
        case 'tcdiff'
            eval([ 'info.tcdiff' apx ' = varargin{j+1};']);
        case 'rsc'
            eval([ 'info.rsc' apx ' = varargin{j+1};']);
        case 'prsc'
            eval([ 'info.prsc' apx ' = varargin{j+1};']);
        case 'rsig'
            eval([ 'info.rsig' apx ' = varargin{j+1};']);
        case 'prsig'
            eval([ 'info.prsig' apx ' = varargin{j+1};']);
        case 'rsc_2nd'
            eval([ 'info.rsc_2nd' apx ' = varargin{j+1};']);
        case 'prsc_2nd'
            eval([ 'info.prsc_2nd' apx ' = varargin{j+1};']);
        case 'rsig_2nd'
            eval([ 'info.rsig_2nd' apx ' = varargin{j+1};']);
        case 'prsig_2nd'
            eval([ 'info.prsig_2nd' apx ' = varargin{j+1};']);
        case 'resvars'
            eval([ 'info.resvars' apx ' = varargin{j+1};']);
        case 'expduration'
            eval([ 'info.expduration' apx ' = varargin{j+1};']);
        case 'sdfs'
            eval([ 'info.sdfs' apx ' = varargin{j+1};']);
        case 'times'
            eval([ 'info.times' apx ' = varargin{j+1};']);
        case 'resdur'
            eval([ 'info.dur' apx ' = varargin{j+1};']);
        case 'ed'
            info.ed = varargin{j+1};
        case 'eX'
            info.eX = varargin{j+1};
        case 'eY'
            info.eY = varargin{j+1};
        case  'phasesel'
            eval([ 'info.phasesel' apx ' = varargin{j+1};']);
        case  'psth'
            eval([ 'info.psth' apx ' = varargin{j+1};']);
        case 'p_anova'
             eval([ 'info.p_anova' apx ' = varargin{j+1};']);
        case 'nrep'
            eval([ 'info.nrep' apx ' = varargin{j+1};']);
        case 'c0rate'
            eval([ 'info.c0rate' apx ' = varargin{j+1};']);
        case 'c0geomn'
            eval([ 'info.c0geomn' apx ' = varargin{j+1};']);
        case 'c0geomn_2nd'
            eval([ 'info.c0geomn_2nd' apx ' = varargin{j+1};']);
        case 'trials_c0'
            eval([ 'info.trials_c0' apx ' = varargin{j+1};']);
        case 'trials_c1'
            eval([ 'info.trials_c1' apx ' = varargin{j+1};']);
    end
    j = j+2;
end

end

