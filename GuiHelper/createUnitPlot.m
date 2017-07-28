function dat = createUnitPlot(exinfo, fctX, fctY, spec, clickFigs, hist_flag)
% dat = createUnitPlot(expInfo, fctX, fctY, spec, fig2plot, hist_flag)
%
%exinfo - the result strucure
%fctX   - the name of the metric to be shown on the x axis
%fctY   - the name of the metric to be shown on the y axis
%spec   - other GUI settings including the stimulus the experiment should
%           test and the ocular condition
%clickFigs - the cell array of strings that name the plots to be shown when
%               the data is clicked, see DataPressed
%hist_flag - flag determining whether the histograms should be shown or not
% 
% 
% @CL
% 

%%%-------------------------------------------------- get the relevant data

% get one experiment per unit and drug condition and spike cluster 
exinfo_res = getUnitComp(spec, exinfo);

% retrieve the metric from the exifno structure according to the axis
% specifications fctX and fctY
val = evalMU(fctX, fctY, exinfo_res);
dat.x = val.x;  dat.y = val.y;
dat.xlab = [val.xlab ' '  spec.stimx ' ' spec.eyex];
dat.ylab = [val.ylab ' ' spec.stimy ' ' spec.eyey];
dat.expInfo = val.exinfo;


clearvars val

% remove nan values
isntnan = find(~isnan(dat.x) & ~isnan(dat.y));
dat.x       = dat.x(isntnan);            
dat.y       = dat.y(isntnan);
dat.is5HT   = [dat.expInfo(isntnan).is5HT];  
dat.expInfo = dat.expInfo(isntnan);  
dat.is5HT   = logical(dat.is5HT);


% write the unit id, x and y values, whether it is cluster 2 and whetehr it
% is a drug experiment with 5HT in the comannd window
disp('id    (x y)   isc2 is5HT \n');
for lll = 1:length(dat.x)
  fprintf('%d   (%1.2f   %1.2f)     %d  %d \n', ...
    dat.expInfo(lll).id, ...
    dat.x(lll), ...
    dat.y(lll), ...
    dat.expInfo(lll).isc2, ...
    dat.expInfo(lll).is5HT);
end


%%%---------------------------------------------------- plot the results
if hist_flag % if the histograms are shown, the main plot is smaller
    pos_hist = [ [0.45 .7 .28 .15];
        [0.78 .18 .12 .45];
        [0.45 .18 .28 .45] ];
else
    pos_hist(3,:)  = [0.45 .18 .5 .7];
end

% dinstinguish the two distributions (5HT and NaCl) via logical indexing
i_5HT = [dat.is5HT]==1;

% specify the range
rx1 = min(dat.x) - abs(min(dat.x)/2);               rx2 = max(dat.x) + abs(max(dat.x)/2);
ry1 = min(dat.y(1, :)) - abs(min(dat.y(1,:))/2);    ry2 = max(dat.y(1, :)) + abs(max(dat.y(1, :))/2);

% specify marker appearance
markerface = zeros(length(dat.x), 3);% (0,0,0) is black
markerface(i_5HT, 1) =  1; % (1,0,0) is red



if hist_flag
    %%%-------------------------------- histogram x
    subplot(3,3, [1 2]);
    plotHist(dat.x, i_5HT);
    set(gca, 'Position', pos_hist(1,:));
    if (rx1 - rx2) ~= 0
        set(gca, 'Xlim', [rx1 rx2])
    end
    %%%--------------------------------- histogram y
    subplot(3,3, [6 9])
    plotHist(dat.y, i_5HT);
    set(gca, 'Position', pos_hist(2, :));
    if (ry1 - ry2) ~= 0
        set(gca, 'Xlim', [ry1 ry2])
    end
    set(gca, 'view',[90 -90]);
    tt = findobj(gca, 'Type', 'text');
    set(tt, 'rotation', -90);
    %%%----------------------------------- main plot
    subplot(3,3, [4 5 7 8])
end


% scatter data in the main plot
% each data point has to be plotted individually to allow the specific
% Callback function that opens the individual plots upon mouse click
for i = 1:length(dat.x)
        scatter(dat.x(i), dat.y(i), 50, ...
            markerAssignment(dat.expInfo(i).param1, dat.expInfo(i).monkey ),...
            'MarkerFaceColor', markerFaceAssignment( dat.expInfo(i) ),...
            'MarkerEdgeColor', markerface(i,:), ...
            'MarkerFaceAlpha', 0.4,...
            'ButtonDownFcn', { @DataPressed, dat.expInfo(i), ...
            clickFigs} );

    hold on;
end


% annotations
xlabel(dat.xlab); ylabel(dat.ylab);

if ~hist_flag
    addTitle(dat);
end

% axes modifications
hold off;
set(gca, 'Position', pos_hist(3,:), 'UserData', dat); % attach the dat strucure to the plot
if (rx1 - rx2) ~= 0
    set(gca, 'xlim',  [rx1 rx2], 'ylim', [ry1 ry2]); % adjust the axes' limits
end


end


