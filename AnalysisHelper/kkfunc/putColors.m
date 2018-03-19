function putColors
% put corresponding colors to the opend figure made by Corinna
% Note that this function does not work for figures which did not derive from Corinna's GUI.
% This function requires Corinna's functions, so the path containing them must be added beforehand by
% 'addpath(genpath(...your path comes here...))'.

% Corinna's function to extract the datapoint
dat = get(gca,'UserData');
len_d = length(dat.x);

% datapoints
points = get(gca, 'Children');
px = nan(1,len_d);
py = nan(1,len_d);
q = 1;
for i = 1:length(points)
    if strcmp(points(i).Type, 'scatter')
        px(q) = [points(i).XData];
        py(q) = [points(i).YData];
        q = q + 1;
    else
        continue
    end
end

% change colors
for i = 1:len_d
  % coresponding color
  col = getCol4Stim(dat.is5HT(i), dat.expInfo(i).param1);

  % search matching datapoint
  c = 100;
  p = 1;
  for k = 1:len_d
    cand = (dat.x(i) - px(k))^2 + (dat.y(i) - py(k))^2;
    if cand < c
        c = cand;
        p = k;
    end
  end
  
  % change colors in the figure
  if dat.expInfo(i).pmodulation(2) < 0.05
    set( points(p), 'MarkerFaceColor', col, 'MarkerEdgeColor', col )
  else
    set( points(p), 'MarkerFaceColor', 'w', 'MarkerEdgeColor', col )
  end
end


