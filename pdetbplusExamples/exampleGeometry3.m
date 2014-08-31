%% Example 3: Create arbitrary region out of arcs and lines
% Create a geometry that is defined by 2 arc boundaries and 1 line boundary. The inner region is
% 'regionA' and the exterior region is 'regionB'. The exterior region is not meshed.

import pdetbplus.*; % import package
% Initial setting for mesh display
showMesh = true; % show mesh in plots

% Arc 1
% Define center of arc
center = pointObject(0.76,-0.5);
% Define start and end points of arc
pt1 = pointObject(0.75,0.5);
pt2 = pointObject(0.85,0.5);
% Define arc
la{1} = arcObject('l1','center',center,'startPoint',pt1,'endPoint',pt2);

% Arc 2
% Define center
center = pointObject(0,0);
% Define start point
startPoint = la{1}.endPoint();
% Define angle of rotation relative to start point
rotationAngle = -pi/8; % negative for clockwise
% Define arc
la{2} = arcObject('l2','center',center,'startPoint',startPoint,'rotationAngle',rotationAngle);

% Line joining end point of arc 2 and start point of arc 1
la{3} = lineObject('l3',la{2}.endPoint(),la{1}.startPoint());

% By setting left and right regions of boundaries, the inner and exterior regions are defined
for k=1:3
    la{k}.leftRegion = 'regionB';
    la{k}.rightRegion = 'regionA';
end

% Define geometry
e3 = geometryObject('arclineShape',la);

% Specify exterior region.
% The exterior region is not meshed.
e3.exteriorRegion = 'regionB';

% Rotate a bit around (0,0)
e3 = e3.rotate(pointObject(0,0),-pi/50);

% Make a copy so e3 can be used elsewhere
tmp = e3;

% Mesh
tmp.initMesh('showMesh',showMesh);
%%
% See help for <matlab:doc('pdetbplus.pointObject') pointObject> <matlab:doc('pdetbplus.arcObject') arcObject> <matlab:doc('pdetbplus.lineObject') lineObject> <matlab:doc('pdetbplus.geometryObject') geometryObject>