classdef arcObject < pdetbplus.boundaryObject
    %arcObject Create arc boundary
    %   arcObject is used for creating geometryObject instances.
    %
    %   arcObject methods:
    %
    %       startPoint      - Get start point of arc as pointObject
    %       endPoint        - Get end point of arc as pointObject
    %       scale           - Scale arc radius and center and return arc
    %       translate       - Translate arc center and return arc
    %       rotate          - Rotate arc with angle around origin and return arc
    %       getXY           - Return XY coordinates as [x,y], given angles
    %       isOn            - Return true if point is on arc, false otherwise
    %
    %   arcObject properties:
    %
    %       startParam      - Value of angle for start point of arc
    %       endParam        - Value of angle for end point of arc
    %       name            - Name of arc
    %       leftRegion      - Name of left region of arc going from start to end points
    %       rightRegion     - Name of right region of arc going from start to end points
    %
    %   Example: create an arc named 'A' with center (1,-1), radius, 5 and startAngle = pi/6 and
    %   endAngle = pi/4.
    %
    %       l1 = arcObject('A','center',pointObject(1,-1),'radius',5,'startAngle',pi/6,'endAngle',pi/4);
    %
    %   See also boundaryObject arcObject piecewiseLineObject parametricLineObject pointObject geometryObject.
    %
    % 
    properties
        %center Center of arc as pointObject
        center;
        %radius Radius of arc
        radius;
    end
    properties (Constant = true,Hidden = true)
        caseNone = 0;
        caseAngleAngle = 1;
        caseStartPointEndPoint = 2;
        caseStartPointRotationAngle = 3;
        caseGradCalculation = 4;
    end
    methods
        % constructor has below name-value formats
        % (name,center,radius,startAngle,endAngle)
        % e.g. ('foo','center',pointObject(0,0),'radius',5,'startAngle',0,pi/2)
        % Note: if endAngle > startAngle, then the arc is anti-clockwise
        % (name,center,startPoint,endPoint)
        % e.g. ('foo','center',pointObject(0,0),'startPoint',pointObject(0,5),'endPoint',pointObject(5,0))
        % (name,center,startPoint,rotationAngle)
        % e.g. ('foo','center',pointObject(0,0),'startPoint',pointObject(0,5),'rotationAngle',pi/2)
        % Note: if rotationAngle > 0, then the arc is anti-clockwise
        % (name, startPoint, endPoint, indexOfPointForGradientCalculation(value of 1 or 2), directionOrthogonalToGradient)
        % e.g. ('foo','startPoint',pointObject(0,5),'endPoint',pointObject(5,0),'indexOfPointForGradientCalculation',1,'directionOrthogonalToGradient',[0 1])
        function obj = arcObject(name,varargin)
            %arcObject Create arc
            %   arcObject(name,varargin) creates arcObject using 4 types of input.
            %       Angle-Angle                         : arcObject(name,'startAngle',startAngleValue,'endAngle',endAngleValue,'center',centerValue,'radius',radiusValue)
            %       Start-End points                    :
            %       arcObject(name,'startPoint',startPointValue,'endPoint',endPointValue,'center',centerValue)
            %       Start Point-Rotation Angle          : arcObject(name,'startPoint',startPointValue,'rotationAngle',rotationAngleValue,'center',centerValue)
            %       Start Point-End Point-Gradient info :
            %       arcObject(name,'startPoint',startPointValue,'endPoint',endPointValue,'indexOfPointForGradientCalculation',indexValue,'directionOrthogonalToGradient',directionValue)
            %
            %   radius refers to the arc radius and center of type pointObject refers to the arc
            %   center.
            %
            %   startPoint and endPoint refer to start and end points of the arc as pointObjects.
            %
            %   startAngle and endAngle are the angles in radians for the start and end points. If
            %   endAngle > startAngle, then the arc is anti-clockwise.
            %
            %   rotationAngle is the difference between endAngle and startAngle in radians. If
            %   rotationAngle > 0, then the arc is anti-clockwise.
            %
            %   In the Start Point-End Point-Gradient info case, the center point and the radius of
            %   the arc are inferred by selecting either the start or end points
            %   (indexOfPointForGradientCalculation) and making the gradient of the arc at that
            %   point orthogonal to the specified direction (directionOrthogonalToGradient).
            %   indexOfPointForGradientCalculation = 1 selects startPoint and
            %   indexOfPointForGradientCalculation = 2 selects endPoint.
            %   directionOrthogonalToGradient is a 2x1 double array.
            %
            %   Examples: create arcs named 'foo'
            %
            %       l1 = arcObject('foo','center',pointObject(0,0),'radius',5,'startAngle',0,pi/2);
            %       l2 = ...
            %       arcObject(('foo','center',pointObject(0,0),'startPoint',pointObject(0,5),'endPoint',pointObject(5,0));
            %       l3 = ...
            %       arcObject('foo','center',pointObject(0,0),'startPoint',pointObject(0,5),'rotationAngle',pi/2);
            %       l4 =
            %       arcObject('foo','startPoint',pointObject(0,5),'endPoint',pointObject(5,0),'indexOfPointForGradientCalculation',1,'directionOrthogonalToGradient',[0 1]);
            %
            %   See also lineObject geometryObject.
            import pdetbplus.arcObject;
            obj.name = name;
            caseType = arcObject.caseNone;
            otherArcData = [];
            for k=1:2:length(varargin)
                a = varargin(k);
                b = varargin(k+1);
                b = b{1};
                if strcmp(a,'center')
                    obj.center = b;
                elseif strcmp(a,'radius')
                    obj.radius = b;
                elseif strcmp(a,'startAngle')
                    obj.startParam = b;
                elseif strcmp(a,'endAngle')
                    obj.endParam = b;
                    caseType = arcObject.caseAngleAngle;
                elseif strcmp(a,'startPoint')
                    otherArcData.startPoint = b;
                elseif strcmp(a,'endPoint')
                    otherArcData.endPoint = b;
                    if caseType ~= arcObject.caseGradCalculation
                        caseType = arcObject.caseStartPointEndPoint;
                    end
                elseif strcmp(a,'rotationAngle')
                    otherArcData.rotationAngle = b;
                    caseType = arcObject.caseStartPointRotationAngle;
                elseif strcmp(a,'indexOfPointForGradientCalculation')
                    otherArcData.indexOfPointForGradientCalculation = b;
                elseif strcmp(a,'directionOrthogonalToGradient')
                    otherArcData.directionOrthogonalToGradient = b;
                    caseType = arcObject.caseGradCalculation;
                end
            end
            switch caseType
                case arcObject.caseAngleAngle
                    % arc corresponding to center, start and end angles
                case arcObject.caseStartPointEndPoint
                    % arc corresponding to center, start and end points
                    p1 = otherArcData.startPoint;
                    obj.radius = norm([p1.x-obj.center.x;p1.y-obj.center.y]);
                    obj.startParam = atan2(p1.y-obj.center.y,p1.x-obj.center.x);
                    p2 = otherArcData.endPoint;
                    obj.endParam = atan2(p2.y-obj.center.y,p2.x-obj.center.x);
                    if obj.startParam <= 0
                        obj.startParam = obj.startParam + 2*pi;
                    end
                    if obj.endParam <= 0
                        obj.endParam = obj.endParam + 2*pi;
                    end
                case arcObject.caseStartPointRotationAngle
                    % arc corresponding to center, start point and rotation angle
                    p1 = otherArcData.startPoint;
                    obj.radius = norm([p1.x-obj.center.x;p1.y-obj.center.y]);
                    obj.startParam = atan2(p1.y-obj.center.y,p1.x-obj.center.x);
                    obj.endParam = obj.startParam + otherArcData.rotationAngle;
                case arcObject.caseGradCalculation
                    % arc corresponding to start and end points, and point index where the gradient
                    % should be orthogonal to specified direction
                    % If pA, pB are the consecutive points, 2 is the point
                    % index and direction, d, is [d1 d2] then the center and radius of the arc
                    % is calculated such that the gradient at pB is made orthogonal to d.
                    % The shortest arc is always returned.
                    p1 = otherArcData.startPoint;
                    p2 = otherArcData.endPoint;
                    [obj.center,obj.radius] = arcObject.getCenterAndRadius(p1,p2,otherArcData.indexOfPointForGradientCalculation,...
                        otherArcData.directionOrthogonalToGradient);
                    obj.startParam = atan2(p1.y-obj.center.y,p1.x-obj.center.x);
                    obj.endParam = atan2(p2.y-obj.center.y,p2.x-obj.center.x);
                    if obj.startParam <= 0
                        obj.startParam = obj.startParam + 2*pi;
                    end
                    if obj.endParam <= 0
                        obj.endParam = obj.endParam + 2*pi;
                    end
                    if (obj.endParam - obj.startParam) > pi
                       obj.startParam = obj.startParam + 2*pi; 
                    elseif (obj.startParam - obj.endParam) > pi
                        obj.endParam = obj.endParam + 2*pi;
                    end
                otherwise
                    error 'unrecognized/incomplete arguments for arcObject instantiation';
            end
        end
        function obj = rotate(self,origin,angle)
            %rotate Rotate arc with angle around origin and return new arc
            %   lineObject = rotate(Origin,Angle) rotates arc around Origin of type pointObject
            %   through Angle in radians and returns new arcObject instance.
            %
            %   Example: rotate arcObject l1 around (1,1) through pi/4 radians and assign new
            %   arcObject instance to l2.
            %
            %       origin = pointObject(1,1);
            %       l2 = l1.rotate(origin,pi/4);
            %
            import pdetbplus.*;
            rotatedCenter = self.center.rotate(origin,angle);
            sP = self.startPoint();
            rotatedStartPoint = sP.rotate(origin,angle);
            obj = arcObject(self.name,'center',rotatedCenter,'startPoint',rotatedStartPoint,'rotationAngle',self.endParam-self.startParam);
            obj.leftRegion = self.leftRegion;
            obj.rightRegion = self.rightRegion;
%             obj.tag = self.tag;
        end
        function obj = translate(self,shiftXYPt)
            %translate Translate center and return new arc
            %   arcObject = translate(shiftXYPt) translates arc center by shiftXYPt of type
            %   pointObject and returns new arcObject instance.
            %
            %   Example: translate arcObject l1's center by (1,1) and assign new
            %   arcObject instance to l2.
            %
            %       t = pointObject(1,1);
            %       l2 = l1.translate(t);
            %            
            obj = self;
            obj.center = obj.center + shiftXYPt;
        end
        function [x,y] = getXY(obj,angles)
            %getXY Return XY coordinates as [x,y] for input angles
            %   [x,y] = getXY(distances) returns (x,y) coordinates corresponding to angles
            %
            %   Example: Get x,y coordinates in arcObject l1 for angles pi/6,pi/4;
            %
            %       angles = [pi/6,pi/4];
            %       [x,y] = l1.getXY(angles);
            %            
            x = obj.radius.*cos(angles);
            y = obj.radius.*sin(angles);
            x = x + obj.center.x;
            y = y + obj.center.y;
        end
        function p = startPoint(obj)
            %startPoint Get start point of arc as pointObject
            %   p = startPoint() returns the start point of arc as pointObject.
            %
            %   Example: get start point of arcObject l1.
            %
            %       p = l1.startPoint();
            %
            p = pdetbplus.pointObject([obj.radius*cos(obj.startParam),obj.radius*sin(obj.startParam)]) + obj.center;
        end
        function p = endPoint(obj)
            %endPoint Get start point of arc as pointObject
            %   p = endPoint() returns the end point of arc as pointObject.
            %
            %   Example: get end point of arcObject l1.
            %
            %       p = l1.endPoint();
            %
            p = pdetbplus.pointObject([obj.radius*cos(obj.endParam), obj.radius*sin(obj.endParam)]) + obj.center;
        end
        function self = scale(self,scaleFactor)
            %scale Scale arc radius and center and return arc
            %   arcObject = scale(factor) returns a copy of current arcObject instance with its dimensions
            %   scaled by factor.
            %
            %   Example: scale radius and center of arcObject, l1 by 25 and assign to l2.
            %
            %       l2 = l1.scale(25);
            %
            self.radius = scaleFactor*self.radius;
            self.center = self.center.scale(scaleFactor);
        end
        function value = isOn(self,p)
            %isOn Return true if point is on arc, false otherwise
            %   isOn(p) returns true if pointObject p is on arc and false otherwise.
            %
            %   Example: find if pointObject p1 is on arcObject l1.
            %
            %       isOn = l1.isOn(p1);
            %
            r = sqrt((p.x-self.center.x)^2 + (p.y-self.center.y)^2);
            if abs(r-self.radius) > eps*self.radius
                value = false;
            else
                pointParam = atan2(p.y-self.center.y,p.x-self.center.x);
                if pointParam <= 0
                    pointParam = pointParam + 2*pi;
                end
                startParam = self.startParam;
                if startParam <= 0
                    startParam = startParam + 2*pi;
                end
                endParam = self.endParam;
                if endParam <= 0
                    endParam = endParam + 2*pi;
                end
                if (pointParam >= startParam) && (pointParam <= endParam)
                    value = true;
                else
                    value = false;
                end
            end
        end
    end
    methods(Static,Hidden=true)
        function [center,radius] = getCenterAndRadius(p1,p2,appliedPointIndex,tangentDirection)
            x1 = p1.x;
            y1 = p1.y;
            x2 = p2.x;
            y2 = p2.y;
            t1 = tangentDirection(1);
            t2 = tangentDirection(2);
            if appliedPointIndex == 1
                sol.a = eval('-(- t2*x1^2 + 2*t1*x1*y1 - 2*t1*x1*y2 + t2*x2^2 + t2*y1^2 - 2*t2*y1*y2 + t2*y2^2)/(2*(t2*x1 - t2*x2 - t1*y1 + t1*y2))');
                sol.b = eval('(t1*x1^2 - 2*t1*x1*x2 + 2*t2*x1*y1 + t1*x2^2 - 2*t2*x2*y1 - t1*y1^2 + t1*y2^2)/(2*t2*x1 - 2*t2*x2 - 2*t1*y1 + 2*t1*y2)');
                sol.r = abs(eval('(((t1^2 + t2^2)*(x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)^2)/(t2*x1 - t2*x2 - t1*y1 + t1*y2)^2)^(1/2)/2'));
            else
                sol.a = eval('(t2*x1^2 - t2*x2^2 - 2*t1*x2*y1 + 2*t1*x2*y2 + t2*y1^2 - 2*t2*y1*y2 + t2*y2^2)/(2*t2*x1 - 2*t2*x2 - 2*t1*y1 + 2*t1*y2)');
                sol.b = eval('-(t1*x1^2 - 2*t1*x1*x2 - 2*t2*x1*y2 + t1*x2^2 + 2*t2*x2*y2 + t1*y1^2 - t1*y2^2)/(2*(t2*x1 - t2*x2 - t1*y1 + t1*y2))');
                sol.r = abs(eval('(((t1^2 + t2^2)*(x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)^2)/(t2*x1 - t2*x2 - t1*y1 + t1*y2)^2)^(1/2)/2'));
            end
            % How the above equations were derived using the
            % Symbolic Math Toolbox
            %             syms x y a b r real;
            %             z = (x-a)^2 + (y-b)^2 - r^2;
            %             z1 = subs(z,[x y],[x1 y1]);
            %             z2 = subs(z,[x y],[x2 y2]);
            %             gradz = gradient(z,[x,y]);
            %             z31 = subs(gradz,[x y],[x1 y1]);
            %             z32 = subs(gradz,[x y],[x2 y2]);
            %             z41 = z31(1)*t1 + z31(2)*t2;
            %             z42 = z32(1)*t1 + z32(2)*t2;
            %             sol1 = solve(z1==0,z2==0,z41==0,a,b,r,'Real',true);
            %             sol2 = solve(z1==0,z2==0,z42==0,a,b,r,'Real',true);
            %             if appliedPointIndex == 1
            %                 sol.a = subs(sol1.a,[x1 y1 x2 y2 t1 t2],[p1.x p1.y p2.x p2.y tangentDirection(1) tangentDirection(2)]);
            %                 sol.b = subs(sol1.b,[x1 y1 x2 y2 t1 t2],[p1.x p1.y p2.x p2.y tangentDirection(1) tangentDirection(2)]);
            %                 sol.r = subs(sol1.r,[x1 y1 x2 y2 t1 t2],[p1.x p1.y p2.x p2.y tangentDirection(1) tangentDirection(2)]);
            %             else
            %                 sol.a = subs(sol2.a,[x1 y1 x2 y2 t1 t2],[p1.x p1.y p2.x p2.y tangentDirection(1) tangentDirection(2)]);
            %                 sol.b = subs(sol2.b,[x1 y1 x2 y2 t1 t2],[p1.x p1.y p2.x p2.y tangentDirection(1) tangentDirection(2)]);
            %                 sol.r = subs(sol2.r,[x1 y1 x2 y2 t1 t2],[p1.x p1.y p2.x p2.y tangentDirection(1) tangentDirection(2)]);
            %             end
            center = pdetbplus.pointObject(sol.a, sol.b);           
            radius = sol.r;            
        end
    end
end