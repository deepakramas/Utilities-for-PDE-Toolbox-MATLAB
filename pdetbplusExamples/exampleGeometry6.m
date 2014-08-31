%% Scale, translate, rotate geometry
% Scale geometry by 1/2, translate it in the X direction by 2.5, rotate it anti-clockwise by pi/4
% radians around point:(2.6,0) and add modified geometry to original geometry.

import pdetbplus.*; % import package
% Initial setting for mesh display
showMesh = true; % show mesh in plots

% Create geometry that will be modified
exampleGeometry5;

% Modify as described above
e6 = e5.scale(0.5).translate(pointObject(2.5,0)).rotate(pointObject(2.6,0),...
    pi/4);
exteriorRegion = e5.exteriorRegion;
e6 = e5 + e6;

% Before calling initMesh(), a single exteriorRegion must have been defined
e6.exteriorRegion = exteriorRegion;

% Make a copy so e6 can be used elsewhere
tmp = e6;

% Mesh
tmp.initMesh('showMesh',showMesh);
%%
% See help for <matlab:doc('pdetbplus.geometryObject') geometryObject>
%%
% <exampleGeometry5.html Example 5: Create geometry by adding other geometries>