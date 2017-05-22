function exinfo = addNumInExp( exinfo )
% loads the spike sorting tables and reads out the number of 5HT
% experiments within the unit recording for each drug experiment
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
% for each filename fname look for the according entry in the sorting file
% if it is a concatenated file, take the average sorting quality


fname = exinfo.fname_drug;
id = floor(exinfo.id); % #unit
monkey = exinfo.monkey; % animal name
drugname = exinfo.drugname; % can be 5HT or NaCl  

if isempty(strfind(fname, 'all.grating'))
    [dn, dt, dt_cum] = getDeltaNexpHelper(fname, id, monkey, sorting, drugname);
else
    % loop through all original experiments, if it is a concatenated file
    load(fname);
    fname_list = getFnames(ex);
    for kk = 1:length(fname_list)
        [dn_list(kk), dt_list(kk), dt_cum_list(kk)] = ...
            getDeltaNexpHelper(fname_list{kk}, id, monkey, sorting, drugname);
    end
    dn.id = min([dn_list.id]);  dn.session = min([dn_list.session]);
    dt.id = min([dt_list.id]);  dt.session = min([dt_list.session]);
    dt_cum.id = min([dt_cum_list.id]);   dt_cum.session = min([dt_cum_list.session]);
    
end

end


function [dn, dt, dt_cum] = getDeltaNexpHelper(fname, id, monkey, sorting, drugname)
% reads out the number of drug experiments dn and the time passed between
% first drug experiment and the current experiment dt in the recording of
% unit #id

% get all indices indicating this unit's experiments
idx_id = findIdxMatchingID(sorting.matfilename, id, monkey);

% use only the drug applied experiments
idx_id = idx_id( ~cellfun(@isempty, strfind(sorting.iontophoresis(idx_id), drugname)) );

% find the index of this exact experiment under observation
idx_f = findIdxMatchingFname(sorting.matfilename, fname, monkey);


% find the session related drug experiments
% use the recording data as indication of a session begin...
idx_session = find( ~cellfun(@isempty, sorting.recordingdate(1:idx_f)) , 1, 'last');
% identify all the rows up to the current experiment...
idx_session = idx_session:idx_f;
% and extract the drug experiments
idx_session = idx_session( ~cellfun(@isempty, strfind(sorting.iontophoresis(idx_session), drugname)) );


   
fprintf('ex of interest: %1.0f \n 1st ex in unit: %1.0f \n 1st ex in session: %1.0f \n\n', ...
    idx_f, idx_id(1), idx_session(1));

%%% compute the number of drug experiments between the session start and
%%% the current experiment
% for this we only need the difference between the indices accounting
dn.id = sum( idx_id < idx_f );
dn.session = sum( idx_session < idx_f );


%%% compute the time between first 5HT application and the current
%%% experiment
[dt.id, dt_cum.id] = getDeltaT(idx_id);
[dt.session, dt_cum.session] = getDeltaT(idx_session);


% use the same, generalized computation to derive the time passed
% between first experiment, identified by idx, and the experiment of interest
% and the cummulative time of drug application
    function [dt, dt_cum] = getDeltaT(idx)
        
        if idx(1) == idx_f
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
            %%% exeriment
            idx = idx(idx<idx_f);
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
% the concatination of both, i.e. the ful filename.

% file directory
if strcmp(monkey , 'ka')
    fdir = ['Z:\data\kaki\' getStringID(id) '\'];
else
    fdir = ['Z:\data\mango\' getStringID(id) '\'];
end

% file name
fname_short = strrep(fname_short, '.mat', '');
[~, i_end] = regexp(fname_short, getStringID(id, monkey));
fname_ful = [fname_short(1:i_end+1) 'c1_sortLH_' fname_short(i_end+2:end) '_' drugname '.mat'];

% concatenate directory and filename
fname = fullfile(fdir, fname_ful);

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

% if isfield(ex, 'Header')
%     if isfield(ex.Header, 'fileID')
%         t = ex.Header.fileID;
%     elseif isfield(ex.Header, 'Sort')
%         t = ex.Header.Sort.DateTime;
%     end
% else
%     disp('search for alternatives');
% end
t_strt = ex.Trials(1).TrialStart;
t_end = ex.Trials(end).TrialEnd;

% t = datestr(t);
% t = datevec(t);

end


function id = getID(fname)

idx = strfind(fname, '_');

id = str2double(fname(idx(1)+1:idx(2)-1));


end