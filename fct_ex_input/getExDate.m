function date = getExDate( ex )
% extracts the date stamp from the header
 

if isfield(ex.Header, 'fileID')
    date = datenum([ex.Header.fileID]);
elseif isfield(ex, 'fileID')
    date = datenum([ex.fileID]);
else
    % use the first time stamp if files were concatenated
    date = datenum([ex.Header.Headers(1).fileID]);
end



