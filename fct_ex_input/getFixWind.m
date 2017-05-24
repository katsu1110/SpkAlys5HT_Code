function [ fwh, fww ] = getFixWind( ex )
%width and height of the fixation window in degrees of visual angle


% degress of visual angle per pixel
dpp = atan(ex.setup.monitorWidth/2/ex.setup.viewingDistance)*180/pi/(ex.setup.screenRect(3)/2);

fwh = dpp * ex.fix.WinH; % height of fixation window (originally in pixel)
fww = dpp * ex.fix.WinW; % width of fixation window (originally in pixel)

    

end

