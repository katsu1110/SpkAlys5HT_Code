function [gridX, gridY] = getGridPosition( ex )
% [gridX, gridY] = getGridPosition( ex )
% 
% returns the electrode position on the surface of the scalp in x and y. 
% 
% 
% @CL 


% The setup position is not always the recording position. 
% Check the spike sorting file for this.
if isfield(ex.setup, 'gridX')
    gridX = ex.setup.gridX;
    gridY = ex.setup.gridY; 
elseif  isfield(ex.setup, 'Left_Hemisphere_gridX')
    gridX = ex.setup.Left_Hemisphere_gridX;
    gridY = ex.setup.Left_Hemisphere_gridY;
end

