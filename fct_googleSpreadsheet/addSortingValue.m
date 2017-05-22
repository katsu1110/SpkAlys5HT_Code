function exinfo = addSortingValue( exinfo )
% loads the spike sorting files and extracts the sorting category and
% assigns it to the corresponding data in the exinfo struct.
%
% Note that concatenated files are traeted individually and are assigned with the mean 
% sorting quality. 
% If there is no sorting quality entry, it is set to 6.
%
% @CL December 07, 2016
% changed @CL February 2, 2017


% load spike sorting files
A = readSpikingTable('Kaki_SpikeSortingFile.csv');
B = readSpikingTable('Mango_SpikeSortingFile.csv');

% retract relevant information
sorting.matfilename = [{B.matfilename}'; {A.matfilename}'];
sorting.isoQc2 = [[B.isoQc2]'; [A.isoQc2]'];
sorting.isoQc1 = [[B.isoQc1]'; [A.isoQc1]']';

% remove empty entries
idx = find(~cellfun(@isempty, sorting.matfilename) & ~cellfun(@(x) strcmp(x, 'NAN'), sorting.matfilename));
sorting.matfilename = sorting.matfilename(idx);
sorting.isoQc2 = sorting.isoQc2(idx);
sorting.isoQc1 = sorting.isoQc1(idx);

% loop through exinfo and assign the spike isolation quality
for i = 1:length(exinfo)
    exinfo(i).spkqual_base = getSpkQ(exinfo(i), exinfo(i).fname, sorting);
    exinfo(i).spkqual_drug = getSpkQ(exinfo(i), exinfo(i).fname_drug, sorting);
end
end


%%
function spkQ = getSpkQ(exinfo, fname, sorting)
% for each filename fname look for the according entry in the sorting file
% if it is a concatenated file, take the average sorting quality

    if isempty(strfind(fname, 'all.grating'))
        spkQ = getSpkQHelper(fname, exinfo.monkey, sorting);
    else
        load(fname);
        fname_list = getFnames(ex);
        for kk = 1:length(fname_list)
            spkQ_list(kk) = getSpkQHelper(fname_list{kk}, exinfo.monkey, sorting);
        end
        spkQ = mean(spkQ_list);
    end 
    
end


function spkQual = getSpkQHelper(fname, monkey, sorting)
% reads out the spike quality saved in the excel file

idx = findIdxMatchingFname(sorting.matfilename, fname, monkey);

if ~isempty(idx)
    if isempty(strfind(fname, 'c1'))
        spkQual = sorting.isoQc2(idx);
    else
        spkQual = sorting.isoQc1(idx);
    end
else
    spkQual = 6;  % when there is no entry, assign 6
end
    
end

