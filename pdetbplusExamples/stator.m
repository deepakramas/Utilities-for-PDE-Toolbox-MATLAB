% create stator
%  
function [assembly] = stator(statorRadius)
import pdetbplus.*;
center = pointObject(0.0,0.0);
if nargin < 1
    statorRadius = 58.40;
end
r1 = statorRadius;
toothHeight = 30.2;
r2 = r1 + toothHeight;
statorThickness = 13.8;
r3 = r2 + statorThickness;
angle1 = 20.2*pi/180;
angle2 = angle1*r1/r2;
windingWidth1 = 6;
windingWidth2 = 2;
windingNotch1 = 2;
windingNotch2 = 10;
assembly = geometryObject('stator',[]);
numParts = 8;
airName = 'air';
for k=1:numParts
    c = [];
    % define boundary points starting from left moving clockwise
    p1 = pointObject([0 r2]);
    p1 = p1.rotate(center,360/(2*numParts)*pi/180);
    p2 = pointObject([0 r3]);
    p2 = p2.rotate(center,360/(2*numParts)*pi/180);
    p3 = p2.rotate(center,-360/numParts*pi/180);
    p4 = p1.rotate(center,-360/numParts*pi/180);
    p6 = pointObject([0 r1]);
    p6 = p6.rotate(center,-angle1/2);
    p7 = pointObject([0 r1]);
    p7 = p7.rotate(center,angle1/2);
    p5 = pointObject([0 r2]);
    p5 = p5.rotate(center,-angle2/2);
    p8 = pointObject([0 r2]);
    p8 = p8.rotate(center,angle2/2);
    
    % 6x2 more points for the winding
    if k == 1 || k == 5
        p9 = p7 + pointObject([0 windingNotch1]);
        p10 = p9 + pointObject([-windingWidth1 0]);
        p11 = p10 + pointObject([0 windingNotch2]);
        p12 = p11 + pointObject([-windingWidth2 0]);
        tmpy = abs(sqrt(r2^2 - p12.x^2));
        p13 = pointObject([p12.x tmpy]);
        p14 = p13 + pointObject([(windingWidth1+windingWidth2) 0]);
        
        p9r = p9.reflect(0,0,0,1);
        p10r = p10.reflect(0,0,0,1);
        p11r = p11.reflect(0,0,0,1);
        p12r = p12.reflect(0,0,0,1);
        p13r = p13.reflect(0,0,0,1);
        p14r = p14.reflect(0,0,0,1);
    end
    % outer
    name = strcat('l1',num2str(k));
    c{end+1} = arcObject(name,'center',center,'startPoint',p2,'endPoint',p3);
    c{end}.leftRegion = 'outer';
    c{end}.rightRegion = 'stator';
    % inner
    if k ~= 1 && k ~= 5
        name = strcat('l2',num2str(k));
        c{end+1} = arcObject(name,'center',center,'startPoint',p4,'endPoint',p5);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        % inner
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p5,p6);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
    
    else
        if k == 1
            windingName = 'winding2';
        else
            windingName = 'winding4';
        end
        name = strcat('l2',num2str(k));
        c{end+1} = arcObject(name,'center',center,'startPoint',p4,'endPoint',p13r);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        % inner
        name = strcat('l3',num2str(k));
        c{end+1} = arcObject(name,'center',center,'startPoint',p13r,'endPoint',p5);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p5,p14r);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p14r,p9r);
        c{end}.leftRegion = windingName;
        c{end}.rightRegion = 'stator';
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p9r,p6);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p14r,p13r);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p13r,p12r);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p12r,p11r);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p11r,p10r);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p10r,p9r);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
    end
    % inner
    name = strcat('l4',num2str(k));
    c{end+1} = arcObject(name,'center',center,'startPoint',p6,'endPoint',p7);
    c{end}.leftRegion = airName;
    c{end}.rightRegion = 'stator';
    % inner
    if k ~= 1 && k ~= 5
        name = strcat('l3',num2str(k));
        c{end+1} = lineObject(name,p7,p8);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        % inner
        name = strcat('l2',num2str(k));
        c{end+1} = arcObject(name,'center',center,'startPoint',p8,'endPoint',p1);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
    else
        if k == 1
            windingName = 'winding1';
        else
            windingName = 'winding3';
        end
        c{end+1} = lineObject(name,p7,p9);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        c{end+1} = lineObject(name,p9,p10);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        c{end+1} = lineObject(name,p10,p11);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        c{end+1} = lineObject(name,p11,p12);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        c{end+1} = lineObject(name,p12,p13);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        c{end+1} = lineObject(name,p13,p14);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = windingName;
        c{end+1} = lineObject(name,p9,p14);
        c{end}.leftRegion = windingName;
        c{end}.rightRegion = 'stator';
        c{end+1} = lineObject(name,p14,p8);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        c{end+1} = arcObject(name,'center',center,'startPoint',p8,'endPoint',p13);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
        c{end+1} = arcObject(name,'center',center,'startPoint',p13,'endPoint',p1);
        c{end}.leftRegion = airName;
        c{end}.rightRegion = 'stator';
    end
    % add an outer domain where mesh terminates
    r4 = r3 + 100;
    pterm = pointObject([0 r4]);
    ptermleft = pterm.rotate(center,360/(2*numParts)*pi/180);
    ptermright = ptermleft.rotate(center,-360/numParts*pi/180);
    c{end+1} = arcObject('outername','center',center,'startPoint',ptermleft,'endPoint',ptermright);
    c{end}.leftRegion = 'nonMeshedSpace';
    c{end}.rightRegion = 'outer';
    part = geometryObject('part',c);
    part = part.rotate(center,2*pi/numParts*(k-1));
    assembly = assembly + part;
end

assembly.exteriorRegion = 'nonMeshedSpace';
end