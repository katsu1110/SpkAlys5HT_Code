function idx = findIdxMatchingFname(fname_spksort, fname, monkey)
% find the index that matches the filename in exinfo with the ones from the
% spike sorting file
% 
%
% fname_spksort is the cell array containing the file names from the
% sorting table
% 
% fname is the file name of the corresponding ex file saved in exinfo.fname /  exinfo.fname_drug 
% 
% monkey is a string specifying 'ma' (mango) or 'ka' (kaki)
%
%
% @CL February 2, 2017



[fnamePref, fnameSuf] = getIDX(fname_spksort, monkey);

idx = find( ~cellfun(@(x) isempty(strfind(fname,x)), fnamePref) & ...
    ~cellfun(@(x) isempty(strfind(fname,x)), fnameSuf));

% sanity check
if length(idx)>2
    fprintf('#%1.0f matching entries found \n', length(idx))
end

% there is few incident where individual experiments were started
% within the same minute. In this case, take the last one. This assumes
% that the previous ones were terminated prematurely
idx = idx(end);

end


function [fnamePref, fnameSuf] = getIDX(fname, monkey)
% splits the name of the recording file in two parts
% this allows to compare the filename of the ex fils with those of the
% trellis (i.e. those saved in the excel file)

for i = 1:length(fname)
    std = regexp(fname{i}, [monkey '_']);
    
    fnamePref{i} = fname{i}(1:7);       % monkey and session number
    
    % try PM as indicator of recording time
    i_strt = regexp(fname{i}, 'PM');
    
    % try AM
    if isempty(i_strt)
        i_strt = regexp(fname{i}, 'AM');
    end
    
    % try .grating
    if isempty(i_strt)
        i_strt = regexp(fname{i}, '.grating');
    end
    
    
    if isempty(i_strt)
        fnameSuf{i} = 'xxx';   % dummy, mostly XPos, YPos or TF experiments
    else
        fnameSuf{i} = fname{i}(i_strt-5:i_strt+2);   % time of recording start
    end
    
    
end
end