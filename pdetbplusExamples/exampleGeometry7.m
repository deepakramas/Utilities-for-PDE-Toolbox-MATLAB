%% Example 7: Create geometry from parametric segments (requires Symbolic Math Toolbox)
% Create geometry named 'paramBased' using 3 parametric boundary segments with inner region called 'regionA' and
% exterior region called 'regionB'. The exterior region is not meshed.
% This functionality depends on a license of the Symbolic Math Toolbox.

import pdetbplus.*; % import package
% Initial setting for mesh display
showMesh = true; % show mesh in plots

% Define parametric curve in terms of parameter r
syms r;
xSym = 2*cos(r);
ySym = 3*sin(r);
% Create segment 1 with parametric curve(2*cos(r),3*sin(r))
le{1} = parametricLineObject('pl4',xSym,ySym,0,pi);
le{1}.leftRegion = 'regionA';
le{1}.rightRegion = 'regionB';

% Modify curve
ySym = ySym + sin(7*r);
% Create segment 2 with modified curve
% Make sure that the second segment starts where the first segment ends
le{2} = parametricLineObject('pl5',xSym,ySym,pi,3/2*pi);
le{2}.leftRegion = 'regionA';
le{2}.rightRegion = 'regionB';

% Create segment 3 - straight line that connects the start and end points of the
% previous segments
le{3} = lineObject('pl6',le{2}.endPoint(),le{1}.startPoint());
le{3}.leftRegion = 'regionA';
le{3}.rightRegion = 'regionB';

% Create geometry
e7 = geometryObject('paramBased',le);

% Specify exterior region.
% The exterior region is not meshed.
e7.exteriorRegion = 'regionB';

% Make a copy so e7 can be used elsewhere
tmp = e7;

% Mesh
tmp.initMesh('showMesh',showMesh);
%%
% See help for <matlab:doc('pdetbplus.parametricLineObject') parametricLineObject> <matlab:doc('pdetbplus.geometryObject') geometryObject>