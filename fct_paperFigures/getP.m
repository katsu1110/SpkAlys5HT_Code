function [p_tt, p_wilc] = getP(ax)
% get p values from the axes title 

s = ax.Title.String{1};
[~,idxptt]  = regexp(s, 'p_{tt}='); idxptt = idxptt+1;
p_tt = str2double( s(idxptt:idxptt+5) );

[~,idxptt]  = regexp(s, 'p_{wilc}='); idxptt = idxptt+1;
p_wilc = str2double( s(idxptt:idxptt+5) );

end
