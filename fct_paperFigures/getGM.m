function gm = getGM(ax)
% get geometric mean from the axes title

[~, idx] = regexp(ax.Title.String, '\mu_{geom}=');
gm(1) = str2double( ax.Title.String{2}(idx{2}+1:end) );
gm(2) = str2double( ax.Title.String{3}(idx{3}+1:end) );

end