% create rotor with angle of rotation = rotateAngle
%  
function assembly = rotor(rotorRadius, airName)
import pdetbplus.*;
if nargin < 1
    rotorRadius = 58.05;
end 
% define point
center = pointObject([0.0 0.0]);
r = rotorRadius;
angle1 = 23.5*pi/180;
% define geometry object
assembly = geometryObject('rotor',[]);
innerArc = cell(0);
for k=1:6
    c = [];
    % define boundary points starting from left moving clockwise
    pomid= pointObject([0 r]);
    poleft = pomid.rotate(center,angle1/2);  
    pileft = poleft + pointObject([0,-16]);
    pimid = pointObject([0 (r-22)]);
    picore = pointObject([0 (r-30)]);
    picoreright = picore.rotate(center,-360/12*pi/180);
    pangleleft = pimid.rotate(center,360/12*pi/180);
    pangleright = pimid.rotate(center,-360/12*pi/180);
    % define outer arc segment
    name = strcat('a1',num2str(k));
    c{end+1} = arcObject(name,'startPoint',pangleleft,'endPoint',pileft,'indexOfPointForGradientCalculation',2,'directionOrthogonalToGradient',[0 1]);
    c{end}.leftRegion = airName;
    c{end}.rightRegion = 'rotor';
    % outer
    name = strcat('l1',num2str(k));
    c{end+1} = lineObject(name,pileft,poleft);
    c{end}.leftRegion = airName;
    c{end}.rightRegion = 'rotor';
    % outer
    name = strcat('a2',num2str(k));
    c{end+1} = arcObject(name,'center',center,'startPoint',poleft,'rotationAngle',-angle1);
    c{end}.leftRegion = airName;
    c{end}.rightRegion = 'rotor';
    % outer
    poright = c{end}.endPoint();     
    piright = poright + pointObject([0,-16]);
    name = strcat('l2',num2str(k));
    c{end+1} = lineObject(name,poright,piright);
    c{end}.leftRegion = airName;
    c{end}.rightRegion = 'rotor';
    % outer
    name = strcat('a3',num2str(k));
    c{end+1} = arcObject(name,'startPoint',piright,'endPoint',pangleright,'indexOfPointForGradientCalculation',1,'directionOrthogonalToGradient',[0 1]);
    c{end}.leftRegion = airName;
    c{end}.rightRegion = 'rotor';
    % inner
    name = strcat('a4',num2str(k));
    innerArc{k} = arcObject(name,'center',center,'startPoint',picoreright,'rotationAngle',360/6*pi/180);
    innerArc{k} = innerArc{k}.rotate(center,2*pi/6*(k-1));
    innerArc{k}.leftRegion = 'core';
    innerArc{k}.rightRegion = 'rotor';
    % create part
    part = geometryObject('rotor',c);
    part = part.rotate(center,2*pi/6*(k-1));
    assembly = assembly + part;
end
innerHole = geometryObject('core',innerArc);
assembly = assembly + innerHole;
assembly.exteriorRegion = airName;
end