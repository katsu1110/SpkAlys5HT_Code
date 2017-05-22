function stats = copyAxes(s, fname, idx)

if nargin == 2; idx = 1; end

h = openfig(fname, 'invisible');

% save the statistics
if nargout ==1
    dat = h.UserData;
    stats = writeStats2File(dat.Stats, fname);
end

% copy the axes
ax = findobj(h, 'type', 'axes');

copyobj(ax(idx).Children, s); delete(h);

end


