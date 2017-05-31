function fname = getFname(ex)
%fname = getFname(ex)
% 
% returns the experiment filename from the ex.Header
% 
% @CL


if isfield(ex.Header, 'Headers')
    fname = ex.Header.Headers(1).fileName;
elseif isfield(ex.Header, 'Sort')
    fname = ex.Header.Sort.exFileName;
else
    fname = ex.Header.fileName;
end