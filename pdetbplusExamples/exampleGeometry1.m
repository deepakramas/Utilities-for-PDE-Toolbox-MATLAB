%% Example 1: Create a circle
% Create a circle named 'circle1' with center = (1,1) and radius = 0.25. The region inside the
% circle is called 'regionA' and the outer region is called 'regionB'. The exterior region is set to
% be the outer region. The exterior region is not meshed.

import pdetbplus.*; % import package
% Initial setting for mesh display
showMesh = true; % show mesh in plots

% Define center
center = pointObject(1,1);
% Define radius
radius = 0.25;
% Create Circle
e1 = geometryObject.createCircle('name','circle1','center',center,'radius',radius,'innerRegion','regionA',...
    'outerRegion','regionB');

% Make a copy so e1 can be used elsewhere
tmp = e1;

% Mesh the circle with Hmax = 2*pi*radius/20 and plot the meshed geometry
tmp.initMesh('showMesh',true,'Hmax',2*pi*radius/20);
%%
% See help for <matlab:doc('pdetbplus.pointObject') pointObject> <matlab:doc('pdetbplus.geometryObject') geometryObject>