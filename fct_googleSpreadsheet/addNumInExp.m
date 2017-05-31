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
% 
% @CL February 2, 2017
% added the temporal distance to the first experiment in the recording
% session


% load spike sorting files
A = readSpikingTable('Kaki_SpikeSortingFile.csv');
B = readSpikingTable('Mango_SpikeSortingFile.csv');
sorting.recordingdate = [{B.VarName1}'; {A.recordingdate}'];
sorting.matfilename = [{B.matfilename}'; {A.matfilename}'];
sorting.iontophoresis = [{B.iontophoresis}'; {A.iontophoresis}'];


% loop throught the result structure and retrieve the information about
% preceeding events
for i = 1:length(exinfo)
    
    fprintf('WORKING ON ROW %1i, file %1.1f \n', i, exinfo(i).id);
    [dn, dt, dt_cum] = ...
        getDeltaNexp(exinfo(i), sorting);
    exinfo(i).dn_id = dn.id;
    exinfo(i).dn_session = dn.session;
    exinfo(i).dt_id = dt.id;
    exinfo(i).dt_session = dt.session;
    exinfo(i).dt_cum_id = dt_cum.id;
    exinfo(i).dt_cum_session = dt_cum.session;
end


end

%%
function [dn, dt, dt_cum] = getDeltaNexp(exinfo, sorting)
% for each filename fname look for the according entry in the sorting file.
% count the number of experiment with 5HT application before this
% experiment within unit (dn_id) and session (dn_session) as well as the
% absolute time to the first 5HT experiment (dt_id, dt_session) and the
% cummulative time of 5HT application (dt_cum_id, dt_cum_session)


fname = exinfo.fname_drug; % the current file
id = floor(exinfo.id); % unit id
monkey = exinfo.monkey; % animal name
drugname = exinfo.drugname; % can be 5HT or NaCl  

if contains(fname, 'all.grating')
    
    % loop through all original experiments, if it is a concatenated file
    load(fname);
    fname_list = getFnames(ex);
    for kk = 1:length(fname_list)
        [dn_list(kk), dt_list(kk), dt_cum_list(kk)] = ...
            getDeltaNexpHelper(fname_list{kk}, id, monkey, sorting, drugname);
    end
     
    % use the minimum distance between experiments, the minimum time, etc.
    % therefore handling the concatenated experiments as one
    dn.id = min([dn_list.id]);  dn.session = min([dn_list.session]);
    dt.id = min([dt_list.id]);  dt.session = min([dt_list.session]);
    dt_cum.id = min([dt_cum_list.id]);   dt_cum.session = min([dt_cum_list.session]);
    
else
    [dn, dt, dt_cum] = getDeltaNexpHelper(fname, id, monkey, sorting, drugname);
end

end


function [dn, dt, dt_cum] = getDeltaNexpHelper(fname, id, monkey, sorting, drugname)
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
[dt.id, dt_cum.id] = getDeltaT(idx_id);
[dt.session, dt_cum.session] = getDeltaT(idx_session);


% use the same, generalized computation to derive the time passed
% between first experiment, identified by idx, and the experiment of interest
% and the cummulative time of drug application
    function [dt, dt_cum] = getDeltaT(idx)
        
        if idx(1) == idx_cex
            dt = 0; % if this is the first experiment
            dt_cum = 0;
        else
            % first drug experiment
            ex1_fname = fulfillfdir(sorting.matfilename{idx(1)}, monkey, ...
                getID(sorting.matfilename{idx(1)}), drugname);
            load( ex1_fname, 'ex' );
            t1 = getTstart(ex);
            
            % current drug experiment
            fname = strrep(fname, 'J:', 'Z:');
            if strcmp(fname(1:2), 'ma');
                fname = fullfile('Z:\data\mango\', getStringID(id, monkey), fname);
            end
            
            load(fname, 'ex'); t2 = getTstart(ex);
            
            % difference in seconds
            dt = (t2 - t1);
            
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

%%
function fname = fulfillfdir(fname_short, monkey, id, drugname)
% this function computes the correct directory and the filename and returns
% the concatination of both, i.e. the full filename.

% file directory
if strcmp(monkey , 'ka')
    fdir = ['Z:\data\kaki\' getStringID(id) '\'];
else
    fdir = ['Z:\data\mango\' getStringID(id) '\'];
end

% filename
fname_short = strrep(fname_short, '.mat', '');
[~, i_end] = regexp(fname_short, getStringID(id, monkey));
fname_ful = [fname_short(1:i_end+1) 'c1_sortLH_' fname_short(i_end+2:end) '_' drugname '.mat'];

% concatenate directory and filename
fname = fullfile(fdir, fname_ful);


% this is hardcoded alrorithms to check for alternative filenames that
% are valid 
if ~exist(fname, 'file')
    fname = strrep(fname, 'XX', 'CO');
end

if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortHN_' fname_short(9:end) '_' drugname '.mat'];
    fname = fullfile(fdir, fname_ful);
end
if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortLH_' fname_short(9:end) '_' drugname 'Ket.mat'];
    fname = fullfile(fdir, fname_ful);
end
if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortLH_' fname_short(9:end) '_' drugname 'SB.mat'];
    fname = fullfile(fdir, fname_ful);
end

if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortLH_' fname_short(9:end) '_' drugname '+GR.mat'];
    fname = fullfile(fdir, fname_ful);
end

if ~exist(fname, 'file')
    fname_ful = [fname_short(1:8) 'c1_sortLH_' fname_short(9:end) '_' drugname 'GR.mat'];
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