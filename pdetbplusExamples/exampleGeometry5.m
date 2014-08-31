%% Example 5: Create geometry by adding other geometries

import pdetbplus.*; % import package
% Initial setting for mesh display
showMesh = true; % show mesh in plots

% Create geometries that will be added
exampleGeometry1;
exampleGeometry2;
exampleGeometry3;
exampleGeometry4;

% The additive geometries should not have intersecting boundaries
e5 = e1 + e2 + e3 + e4;

% Set exterior region to be the same as that for the polygon i.e. regionA
% Before calling initMesh(), a single exteriorRegion must have been defined
e5.exteriorRegion = e4.exteriorRegion;

% Make a copy so e5 can be used elsewhere
tmp = e5;

% Mesh
tmp.initMesh('showMesh',showMesh);
%%
% See help for <matlab:doc('pdetbplus.geometryObject') geometryObject>
%%
% <exampleGeometry1.html Example 1: Create a circle>
%%
% <exampleGeometry2.html Example 2: Create a square>
%%
% <exampleGeometry3.html Example 3: Create arbitrary region out of arcs and lines>
%%
% <exampleGeometry4.html Example 4: Create a polygon>