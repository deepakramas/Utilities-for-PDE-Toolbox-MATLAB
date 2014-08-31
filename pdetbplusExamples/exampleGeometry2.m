%% Example 2: Create a square
% Create a square named 'square2' with coordinates with first edge as (0.25,0.5)->(0.5,0.65). The
% square is completed by moving in the clockwise direction till the first point i.e. (0.25,0.5) is
% met. The region inside the square is called 'regionA' and the outer region is called 'regionB'.
% The outer region is set to be the exterior region by default. The exterior region is not meshed.

import pdetbplus.*; % import package
% Initial setting for mesh display
showMesh = true; % show mesh in plots

% Define starting edge
pt1 = pointObject(0.25,0.5);
pt2 = pointObject(0.5,0.65);
% Create geometry
e2 = geometryObject.createSquare('name','square2','edgeStartPoint',pt1,...
    'edgeEndPoint',pt2,'innerRegion','regionA','outerRegion','regionB','clockwise',true);

% Make a copy so e2 can be used elsewhere
tmp = e2;

% Mesh
tmp.initMesh('showMesh',showMesh);
%%
% See help for <matlab:doc('pdetbplus.pointObject') pointObject> <matlab:doc('pdetbplus.geometryObject') geometryObject>