%% Example 1: How to create a circle
import pdetbplus.*; % import package for accessing the geometryObject and boundary classes
% Initial setting for mesh display
showMesh = true; % show mesh in plots
% Define center : (1,1) as a pointObject type
center = pointObject(1,1);
% Define radius of 0.25
radius = 0.25;
% Specify region inside circle as 'regionA'
innerRegion = 'regionA';
% Specify region outside circle as 'regionB'
exteriorRegion = 'regionB';
% Create Circle called 'circle1'
tmp = geometryObject.createCircle('name','circle1','center',center,'radius',radius,'innerRegion',innerRegion,...
    'outerRegion',exteriorRegion);
% Mesh the circle with Hmax = 2*pi*radius/20 and plot the meshed geometry
tmp.initMesh('showMesh',true,'Hmax',2*pi*radius/20);