function exinfo = addNumInExp( exinfo )
% exinfo = addNumInExp( exinfo )
%
% loads the spike sorting tables and reads out the number of drug
% experiments within the unit recording for each drug experiment. if NaCl
% was applied in the experiment, it says how many NaCl experiments were
% performed previously and vice versa for experiments with 5HT application.
% (go to the subfunction getDeltaNexpHelper if you want to change that).
%
%
%
%
% Following fields are added:
%   dn_id           - the number of preceeding drug experiments before the
%                       current one during the recording of this unit
%   dn_session      -the number of drug experiments before the current one
%                       during the recording session (including other
%                       unit recordings)
%   dt_id           - the absolut time between this and the beginning of
%                       the first drug experiment in the recording of this
%                       unit before the current experiment
%   dt_session      - the absolut time between this and the beginning of
%                       the first drug experiment in this recording session
%   dt_cum_id       - the cummulated time of drug application in the
%                       recording of this unit before the given experiment
%   dt_cum_session  - the cummulated time of drug application before the
%                       given experiment in the recording session
%
%
% @CL February 2, 2017
% added the temporal distance to the first experiment in the recording
% session
%
% @CL June 27, 2017
% added the number of recovery experiments and the time between the last
% drug application and the baseline experiment
%     t_drugapp     - duration of the drug application
%     n_recov       - number of experiments until recovered
%     t_recov       - time until recovery
%
%
%
%
% @CL 12/07/2017
% added the filename of the preceeding drug experiment to exinfo
%
%
%
%

% load spike sorting files
A = readSpikingTable('Kaki_SpikeSortingFile.csv');
B = readSpikingTable('Mango_SpikeSortingFile.csv');
sorting.recordingdate = [{B.VarName1}'; {A.recordingdate}'];
sorting.matfilename = [{B.matfilename}'; {A.matfilename}'];
sorting.iontophoresis = [{B.iontophoresis}'; {A.iontophoresis}'];
sorting.isoQc1 = [[B.isoQc1]'; [A.isoQc1]']';


% loop throught the result structure and retrieve the information about
% preceeding events
errorrows = [];
parfor i = 1:length(exinfo)
    try
        fprintf('WORKING ON ROW %1i, file %1.1f \n', i, exinfo(i).id);
        exinfo(i) = helpdeltan(exinfo(i), sorting);
    catch
        errorrows = [errorrows, i];
    end
   
end


fprintf('problems with \n %1d rows \n', errorrows);
save('errorrows_addnum.mat', 'errorrows');
end


function exinfo = helpdeltan(exinfo, sorting)

[dn, dt, dt_cum, t_app, recov_baseline, recov_recovery] = ...
    getDeltaNexp(exinfo, sorting);
exinfo.dn_id = dn.id;
exinfo.dn_session = dn.session;
exinfo.dt_id = dt.id;
exinfo.dt_session = dt.session;
exinfo.dt_cum_id = dt_cum.id;
exinfo.dt_cum_session = dt_cum.session;

% information about the recovery of the baseline experiment
exinfo.t_drugapp = t_app;
exinfo.n_recov_base = recov_baseline.n;
exinfo.t_recov_base = recov_baseline.t;
exinfo.recov_base_fname = recov_baseline.prev_drug_fname; % previous drug experiment

% information about the following recovery
exinfo.t_recov = recov_recovery.t;
exinfo.recov_fname = recov_recovery.fname;
exinfo.recov_p = recov_recovery.p;

end

%%
function [dn, dt, dt_cum, t_app, recov_baseline, recov_recovery] = getDeltaNexp(exinfo, sorting)
% for each filename fname look for the according entry in the sorting file.
% count the number of experiment with 5HT application before this
% experiment within unit (dn_id) and session (dn_session) as well as the
% absolute time to the first 5HT experiment (dt_id, dt_session) and the
% cummulative time of 5HT application (dt_cum_id, dt_cum_session)


fname = exinfo.fname_drug; % the current drug experiment filename
id = floor(exinfo.id); % unit id
monkey = exinfo.monkey; % animal name
drugname = exinfo.drugname; % can be 5HT or NaCl

if contains(fname, 'all.grating')
    
    % loop through all original experiments, if it is a concatenated file
    load(fname);
    fname_list = getFnames(ex);
    
    for kk = 1:length(fname_list)
        [dn_list(kk), dt_list(kk), dt_cum_list(kk), t_app(kk)] = ...
            getDeltaNexpHelper(fname_list{kk}, id, monkey, sorting, drugname);
    end
    
    % use the minimum distance between experiments, the minimum time, etc.
    % therefore handling the concatenated experiments as one
    dn.id = min([dn_list.id]);  dn.session = min([dn_list.session]);
    dt.id = min([dt_list.id]);  dt.session = min([dt_list.session]);
    dt_cum.id = min([dt_cum_list.id]);   dt_cum.session = min([dt_cum_list.session]);
    
    % use the summed application duration
    t_app = sum(t_app);
    
    
else
    [dn, dt, dt_cum, t_app] = getDeltaNexpHelper(fname, id, monkey, sorting, drugname);
    
end


%%% -- obtain information regarding the recovery before the baseline experiment

fname = exinfo.fname; % the current baseline experiment filename

if contains(fname, 'all.grating')
    
    % loop through all original experiments, if it is a concatenated file
    load(fname);
    fname_list = getFnames(ex);
    
    for kk = 1:length(fname_list)
        recov_baseline_list(kk) = getRecoveryDuration(fname_list{kk}, id, monkey, sorting, drugname);
        recov_recov_list(kk) = getRecovFname(fname_list{kk}, id, monkey, sorting, exinfo);
        
    end
    
    % use the summed application duration and the minimum distance for the
    % recovery
    recov_baseline.n = min([recov_baseline_list.n]);
    recov_baseline.t = min([recov_baseline_list.t]);
    
    [~, mini] = min([recov_baseline_list.n]);
    recov_baseline.prev_drug_fname = recov_baseline_list(mini).prev_drug_fname;
    
    [~,mini] = min([recov_recov_list.t]);
    recov_recovery.t = recov_recov_list(mini).t;
    recov_recovery.fname = recov_recov_list(mini).fname;
    recov_recovery.p = recov_recov_list(mini).p;
    
    
else
    recov_baseline = getRecoveryDuration(fname, id, monkey, sorting, drugname);
    recov_recovery = getRecovFname(fname,id, monkey, sorting, exinfo);
    
end


end

%%
function recov = getRecoveryDuration(fname, id, monkey, sorting, drugname)
% reads out the number of recovery experiments and the recovery duration
% before the current baseline experiment

% get all indices indicating this unit's experiments in the sorting file
idx_id = findIdxMatchingID(sorting.matfilename, id, monkey);

% find the index of the current experiment under observation
idx_cex = findIdxMatchingFname(sorting.matfilename, fname, monkey);

% only consider experiments before this - reduce idx_id
idx_id = idx_id(1):idx_cex;

% within these, find the last drug exepriment
idx_drug = idx_id( contains(sorting.iontophoresis(idx_id), drugname));
% alternatively, look for the last iontophoretic application
% idx_drug = idx_id( contains(sorting.iontophoresis(idx_id), '5HT') || ...
%           contains(sorting.iontophoresis(idx_id), 'NaCl') );


if isempty(idx_drug)
    recov.n = -1;
    recov.t = -1;
    recov.prev_drug_fname = [];
else
    idx_prev_drug = idx_drug(end);
    
    recov.n = idx_cex-idx_prev_drug-1;
    
    fname = strrep(fname, 'J:', 'Z:');
    if strcmp(fname(1:2), 'ma')
        fname = fullfile('Z:\data\mango\', getStringID(id, monkey), fname);
    end
    exb = load(fname, 'ex'); exb = exb.ex; % current baseline experiment
    
    [t_strt_cex, ~]=getTstart(exb);
    
    
    prev_drug_fname = fulfillfdir(sorting.matfilename{idx_prev_drug}, monkey, ...
        getID(sorting.matfilename{idx_prev_drug}), drugname);
    exd = load(prev_drug_fname, 'ex');    exd = exd.ex;
    
    [~, t_end_drugex]=getTstart(exd);
    
    recov.t = t_strt_cex-t_end_drugex;
    pause(3)
    
    recov.prev_drug_fname = prev_drug_fname;
    
end

end


%% get the filename of the following recovery experiment

function recov = getRecovFname(fname, id, monkey, sorting, exinfo)


recov.fname = [];
recov.t = -1;
recov.p = 1;

fname_drug = exinfo.fname_drug;

if contains(fname_drug, 'all.grating')
    load(fname_drug);
    fname_list = sort(getFnames(ex));
    
    fname_drug = fname_list{1};
end

% the stimulus dimension tested in this experiment
stimdim = getStimDim(fname);

% indices of this unit's experiments in the sorting file
idx_id = findIdxMatchingID(sorting.matfilename, id, monkey);

% the current experiment under observation
idx_cex = findIdxMatchingFname(sorting.matfilename, fname, monkey);


if idx_cex<idx_id(end)
    
    % potential recovery experiments
    idx_recov_list = idx_cex+1:idx_id(end);
    
    % if there is a drug experiment in the present recovery list, short
    % the list up to this experiment
    idx_drug = idx_recov_list( contains(sorting.iontophoresis(idx_recov_list), '5HT') | ...
        contains(sorting.iontophoresis(idx_recov_list), 'NaCl') );
    
    
    if length(idx_drug)>1
        idx_recov_list = idx_drug(1):idx_drug(2);
    elseif ~isempty(idx_drug) 
        idx_recov_list = idx_drug(1):idx_id(end); 
    end
    
    % only control experiments with the same stimulus conditions are considered
    if strcmp(stimdim, 'RC')
        idx_recov = idx_recov_list( contains(sorting.matfilename(idx_recov_list), stimdim) & ...
            contains(sorting.iontophoresis(idx_recov_list), 'control') & ...
            sorting.isoQc1(idx_recov_list)'<4.5 );
    else
        idx_recov = idx_recov_list( contains(sorting.matfilename(idx_recov_list), stimdim) & ...
            ~contains(sorting.matfilename(idx_recov_list), 'RC') & ...
            contains(sorting.iontophoresis(idx_recov_list), 'control') & ...
            sorting.isoQc1(idx_recov_list)'<4.5 );
    end
    
    
    if ~isempty(idx_recov)
        
        % exception for mangos files
        fname = strrep(fname, 'J:', 'Z:');
        if strcmp(fname(1:2), 'ma')
            fname = strrep('Z:\data\mango\', getStringID(id, monkey), fname);
        end
        
        % current baseline experiment
        ex_base = loadCluster(fname, 'ocul', exinfo.ocul, 'loadlfp', false);
        
        % loope through the recovery data and find the first that is not
        % significantly different from the baseline
        for kk = 1:length(idx_recov)
            recov_fname_temp{kk} = fulfillfdir(sorting.matfilename{idx_recov(kk)}, monkey, id);
        end
        
        try
            recov_fname_temp = recov_fname_temp(~contains(recov_fname_temp, 'adapt'));
            if ~isempty(recov_fname_temp)
                for kk = 1:length(recov_fname_temp)
                    
                    % next potential recovery experiment
                    ex_recov = loadCluster(recov_fname_temp{kk}, ...
                        'ocul', exinfo.ocul, 'loadlfp', false);
                    
                    p_recov = pModulation(ex_recov, ex_base, exinfo);
                    if p_recov<0.05
                        break;
                    end
                end
                
                % for the recovery duration, we have to load the drug experiment
                ex_drug = load(fname_drug); ex_drug = ex_drug.ex;
                ex_recov = load(recov_fname_temp{kk}); ex_recov = ex_recov.ex;
                
                [t_strt_rec, ~]=getTstart(ex_recov);
                [~, t_end_cex]=getTstart(ex_drug);
                recov.t = t_strt_rec-t_end_cex;
                recov.fname = recov_fname_temp;
                recov.p = p_recov(2);
                recov
                
            end
        catch
           disp('there is a problem'); 
            
        end
        
    end
    
end




end



%%
function [dn, dt, dt_cum, t_app] = getDeltaNexpHelper(fname, id, monkey, sorting, drugname)
% reads out the number of drug experiments dn and the time passed between
% first drug experiment and the current experiment dt in the recording of
% unit #id


% get all indices indicating this unit's experiments
idx_id = findIdxMatchingID(sorting.matfilename, id, monkey);

% use only the drug applied experiments
idx_id = idx_id( contains(sorting.iontophoresis(idx_id), drugname) );

% find the index of the current experiment under observation
idx_cex = findIdxMatchingFname(sorting.matfilename, fname, monkey);


% find the session related drug experiments
% (sessions means recordings on the same date)
% use the recording data as indication of a session begin...
idx_session = find( ~cellfun(@isempty, sorting.recordingdate(1:idx_cex)) , 1, 'last');
% identify all the rows up to the current experiment...
idx_session = idx_session:idx_cex;
% and extract the drug experiments
idx_session = idx_session( contains(sorting.iontophoresis(idx_session), drugname) );



fprintf('ex of interest: %1.0f \n 1st ex in unit: %1.0f \n 1st ex in session: %1.0f \n\n', ...
    idx_cex, idx_id(1), idx_session(1));

%%% compute the number of drug experiments between the session start and
%%% the current experiment
% for this we only need the difference between the indices
dn.id = sum( idx_id < idx_cex );
dn.session = sum( idx_session < idx_cex );


%%% compute the time between first 5HT application and the current
%%% experiment
% current drug experiment
fname = strrep(fname, 'J:', 'Z:');
if strcmp(fname(1:2), 'ma')
    fname = fullfile('Z:\data\mango\', getStringID(id, monkey), fname);
end
load(fname, 'ex');

[t_strt_cex, t_end_cex] = getTstart(ex); % start of the current experiment
t_app = t_end_cex-t_strt_cex;

[dt.id, dt_cum.id] = getDeltaT(idx_id, t_strt_cex);
[dt.session, dt_cum.session] = getDeltaT(idx_session, t_strt_cex);


% use the same, generalized computation to derive the time passed
% between first experiment, identified by idx, and the experiment of interest
% and the cummulative time of drug application
    function [dt, dt_cum] = getDeltaT(idx, t_strt_cex)
        
        if idx(1) == idx_cex
            dt = 0; % if this is the first experiment
            dt_cum = 0;
        else
            % first drug experiment
            ex1_fname = fulfillfdir(sorting.matfilename{idx(1)}, monkey, ...
                getID(sorting.matfilename{idx(1)}), drugname);
            load( ex1_fname, 'ex' );
            t1 = getTstart(ex); % start of the first experiment
            
            % difference in seconds between the last and the
            dt = (t_strt_cex - t1);
            
            %%% the cummulative time of drug application before the start of this
            %%% experiment
            idx = idx(idx<idx_cex);
            dt_cum = 0;
            for kk =1:length(idx)
                ex_fname = fulfillfdir(sorting.matfilename{idx(kk)}, monkey,...
                    getID(sorting.matfilename{idx(kk)}), drugname);
                load( ex_fname, 'ex' );
                [t_strt, t_end]  = getTstart(ex);
                
                dt_cum = dt_cum + (t_end - t_strt);
            end
        end
    end

end




%% subfunctions
function fname = fulfillfdir(fname_short, monkey, id, drugname)
% this function computes the correct directory and the filename and returns
% the concatination of both, i.e. the full filename.


if nargin == 4
    drugsuffix = ['_' drugname];
else
    drugsuffix = '';
end


% file directory
if strcmp(monkey , 'ka')
    fdir = ['Z:\data\kaki\' getStringID(id) '\'];
else
    fdir = ['Z:\data\mango\' getStringID(id) '\'];
end

% filename
fname_short = strrep(fname_short, '.mat', '');
[~, i_end] = regexp(fname_short, getStringID(id, monkey));
fname_ful = [fname_short(1:i_end+1) 'c1_sortLH_' fname_short(i_end+2:end) drugsuffix '.mat'];


% concatenate directory and filename
fname = fullfile(fdir, fname_ful);


% this is hardcoded alrorithms to check for alternative filenames that
% are valid
if ~exist(fname, 'file')
    fname = strrep(fname, 'XX', 'CO');
end

if ~exist(fname, 'file')
    fname_ful = strrep(fname_ful, '.mat', '_adapt.mat');
    fname = fullfile(fdir, fname_ful);
end


if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortHN_' fname_short(9:end) drugsuffix '.mat'];
    fname = fullfile(fdir, fname_ful);
end

if ~exist(fname, 'file')
    fname_ful = strrep(fname_ful, '.mat', '_adapt.mat');
    fname = fullfile(fdir, fname_ful);
end


if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortLH_' fname_short(9:end) drugsuffix 'Ket.mat'];
    fname = fullfile(fdir, fname_ful);
end
if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortLH_' fname_short(9:end) drugsuffix 'SB.mat'];
    fname = fullfile(fdir, fname_ful);
end

if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortLH_' fname_short(9:end) drugsuffix '+GR.mat'];
    fname = fullfile(fdir, fname_ful);
end

if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortLH_' fname_short(9:end) drugsuffix 'GR.mat'];
    fname = fullfile(fdir, fname_ful);
end


if ~exist(fname, 'file')
    fname
    error('this file does not exist as it is');
end


end


%%
function [t_strt, t_end] = getTstart(ex)
% returns the time stamp of the start of the first trial and the end of the
% last trial

t_strt = ex.Trials(1).TrialStart;
t_end = ex.Trials(end).TrialEnd;

end





function id = getID(exname)
% retracts the unit ID from the ex filename
% exname is the identifier from the google spreadsheet
% example: ma_0322_16.05.grating.SF.mat

idx = strfind(exname, '_');
id = str2double(exname(idx(1)+1:idx(2)-1));


end



function stimdim = getStimDim(fname)
% returns the tested stimulus dimension


if contains(fname, 'OR') && ~contains(fname, 'RC')
    stimdim  = 'OR';
elseif contains(fname, 'CO') && ~contains(fname, 'RC')
    stimdim = 'CO';
elseif contains(fname, 'SF') && ~contains(fname, 'RC')
    stimdim = 'SF';
    
elseif contains(fname, 'SZ') && ~contains(fname, 'RC')
    stimdim = 'SZ';
else
    stimdim = 'RC';
end

end





