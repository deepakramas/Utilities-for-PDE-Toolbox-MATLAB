%% Example 8: Create geometry from image (requires Image Processing Toolbox)
% Read image, identify relevant boundaries and create PDE geometry with inner regions - 'regionA'
% and 'regionB' and exterior region - 'regionC'. The exterior region is not meshed. The geometry is
% further scaled and translated.

import pdetbplus.*; % import package
% Initial setting for mesh display
showMesh = true; % show mesh in plots

% Read image and plot
im = imread('testImage.png');
image(im); snapnow;

% Create geometries from recognized boundaries of image
[imgs,isClockwise] = geometryObject.createGeometriesFromImage('image','testImage.png','minimumPoints',15);

% Notes:
% * Not all geometries created this way are interesting or what we expect. 
% We need to manually verify by plotting the geometries individually.
% Here we combine the ones with indices 3,4,..end that *apriori* have
% been determined to be the ones we are interested in (by calling
% the plot() method per geometry)
% * Set interior and exterior regions for geometry 3, we want regionC to be exterior and regionB to
% be interior
e8 = imgs{3}.setInteriorExteriorRegions('regionB','regionC',isClockwise);
for k=4:length(imgs)
    imgs{k} = imgs{k}.setInteriorExteriorRegions('regionA','regionB',isClockwise);
    e8 = e8 + imgs{k};
end

% Specify exterior region.
% The exterior region is not meshed.
e8.exteriorRegion = 'regionC';

% Get limits of geometry to transform the coordinates
[xmin,ymin,xmax,ymax] = e8.getLimitsXY();
% Get rid of XY offsets and scale model to be of size ~ 1
e8 = e8.translate(pointObject(0-xmin,0-ymin)).scale(1/max(abs(xmax-xmin),abs(ymax-ymin)));

% Make a copy so e8 can be used elsewhere
tmp = e8;

% Mesh
tmp.initMesh('showMesh',showMesh); snapnow;
%%
% See help for <matlab:doc('pdetbplus.geometryObject') geometryObject>