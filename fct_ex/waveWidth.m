function wdt = waveWidth( exinfo, w0, w1, p_flag )



[ wdt(1), mnw0, sdw0, tstrt0] = waveWidth_helper( w0 );
[ wdt(2), mnw1, sdw1, tstrt1 ] = waveWidth_helper( w1 );

exinfo.wdt = wdt;
if p_flag
    wavePlot(mnw0, mnw1, sdw0, sdw1, tstrt0, tstrt1, exinfo);
end

end


function [ wdt, mnwave, sdwave, tstrt ] = waveWidth_helper( Waves )
%WAVEWIDTH returns the duration between lowest and following highest point
%in the mean spike wave form.


lx = size(Waves, 2);
timeX_oq = 1:1700/lx:1700;
timeX = 1:5:1700;


Waves(isinf(Waves)) = nan;

if all(isnan(Waves))
    
    wdt = 0;
    tstrt = 0;
    mnwave = zeros(length(timeX),1);
    sdwave = zeros(length(timeX),1);
else
    
    % mean and standard deviation
    mnwave = interp1(timeX_oq, nanmean(Waves), timeX, 'spline');
    sdwave = interp1(timeX_oq, nanstd(Waves), timeX, 'spline');
    
    
    % wave width - distance between peak and trough
    [~, idxmin] = findpeaks(-mnwave,'MinPeakHeight', 0.025);
    idxmin = idxmin(1);
    
    [~, idxmax] = max(mnwave);   % find first maximum after min
    
    if idxmin > idxmax
        temp = idxmin;
        idxmin = idxmax;
        idxmax = temp;
    end
    
    tstrt =timeX(idxmin);
    wdt = timeX(idxmax) - timeX(idxmin);
    % mnwave = mnwave / abs(mnwave(idxmax) - mnwave(idxmin));
end
end





function wavePlot(mnw0, mnw1, sdw0, sdw1, tstrt0, tstrt1, exinfo)
%%% plot

c = getCol(exinfo);
timeX = 1:5:1700;

h = figure('Name', exinfo.figname);
plot(timeX, mnw0, 'b', 'LineWidth', 3); hold on
plot(timeX, mnw1, c, 'LineWidth', 3);
plot(timeX, mnw0-sdw0, timeX, mnw0+sdw0, 'b', 'LineStyle', ':');
plot(timeX, mnw1-sdw1, timeX, mnw1+sdw1, c, 'LineStyle', ':');
plot([tstrt0 tstrt0], get(gca, 'YLim'), 'b--');
plot([tstrt1 tstrt1], get(gca, 'YLim'), '--', 'Color', c);


legend('baseline', exinfo.drugname);

% ylim([-0.4,0.4]);

ylabel('amplitude');
xlabel('time in \mus');
title(sprintf('mean wave form (wdt = %1.2f / %1.2f )', exinfo.wdt));

savefig(h, exinfo.fig_waveform);
close(h);
end
