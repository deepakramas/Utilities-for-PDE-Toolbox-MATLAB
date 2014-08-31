classdef lineObject < pdetbplus.boundaryObject
    %lineObject Create straight line boundary
    %   lineObject is used for creating geometryObject instances.
    %
    %   lineObject methods:
    %
    %       startPoint      - Get start point of line as pointObject
    %       endPoint        - Get end point of line as pointObject
    %       scale           - Scale line dimensions and return line
    %       translate       - Translate line and return line
    %       rotate          - Rotate line with angle around origin and return line
    %       getXY           - Return XY coordinates as [x,y], given param value along line
    %       isOn            - Return true if point is on line, false otherwise
    %
    %   lineObject properties:
    %
    %       startParam      - Value of parameter for start point of line
    %       endParam        - Value of parameter for end point of line
    %       name            - Name of line
    %       leftRegion      - Name of left region of line going from start to end points
    %       rightRegion     - Name of right region of line going from start to end points
    %
    %   Example: create a line named 'A' between (1,1) and (5,3.3) and assign to l1.
    %
    %       p1 = pointObject(1,1); p2 = pointObject(5,3.3);
    %       l1 = lineObject('A',p1,p2);
    %
    %   See also boundaryObject arcObject piecewiseLineObject parametricLineObject pointObject geometryObject.
    %
    % 
    properties
        %p1 Start point as type pointObject
        p1;
        %p2 End point as type pointObject
        p2;
    end
    methods
        % constructor has input types (name, pointObject, pointObject)
        function obj = lineObject(name,point1,point2)
            %lineObject Create lineObject
            %   lineObject(Name,p1,p2) creates lineObject instance named Name given start and end pointObjects p1
            %   and p2 respectively.
            %
            %   Example: create a line named 'A' between (1,1) and (5,3.3) and assign to l1.
            %
            %       p1 = pointObject(1,1); p2 = pointObject(5,3.3);
            %       l1 = lineObject('A',p1,p2);
            %
            obj.p1 = point1;
            obj.p2 = point2;
            obj.startParam = 0;
            obj.endParam = norm([point1.x-point2.x;point1.y-point2.y]);
            obj.name = name;
        end 
        function obj = rotate(self,origin,angle)
            %rotate Rotate line with angle around origin and return line
            %   lineObject = rotate(Origin,Angle) rotates line around Origin of type pointObject
            %   through Angle in radians and returns new lineObject instance.
            %
            %   Example: rotate lineObject l1 around (1,1) through pi/4 radians and assign new
            %   lineObject instance to l2.
            %
            %       origin = pointObject(1,1);
            %       l2 = l1.rotate(origin,pi/4);
            %
            import pdetbplus.*;
            pt1 = self.p1.rotate(origin,angle);
            pt2 = self.p2.rotate(origin,angle);
            obj = lineObject(self.name,pt1,pt2);
            obj.leftRegion = self.leftRegion;
            obj.rightRegion = self.rightRegion;
        end
        function obj = translate(self,shiftXYPt)
            %translate Translate line and return line
            %   lineObject = translate(shiftXYPt) translates line around shiftXYPt of type
            %   pointObject and returns new lineObject instance.
            %
            %   Example: translate lineObject l1 through (1,1) and assign new
            %   lineObject instance to l2.
            %
            %       t = pointObject(1,1);
            %       l2 = l1.translate(t);
            %
            obj = self;
            obj.p1 = obj.p1 + shiftXYPt;
            obj.p2 = obj.p2 + shiftXYPt;
        end
        function [x,y] = getXY(obj,distances)
            %getXY Return XY coordinates as [x,y], given param value along line
            %   [x,y] = getXY(distances) returns (x,y) coordinates for distances from start to end
            %   points of lineObject instance.
            %
            %   Example: Get x,y coordinates for distances 1,3 from start to end points of lineObject l1;
            %
            %       distances = [1,3];
            %       [x,y] = l1.getXY(distances);
            %
            x = (obj.p2.x - obj.p1.x)*distances/obj.endParam + obj.p1.x;
            y = (obj.p2.y - obj.p1.y)*distances/obj.endParam + obj.p1.y;
            if ~isempty(find(distances == obj.startParam,1))
                x(distances == obj.startParam) = obj.p1.x;
                y(distances == obj.startParam) = obj.p1.y;
            elseif ~isempty(find(distances == obj.endParam,1))
                x(distances == obj.endParam) = obj.p2.x;
                y(distances == obj.endParam) = obj.p2.y;
            end                
        end
        function p = startPoint(obj)
            %startPoint Get start point of line as pointObject
            %   p = startPoint() returns the start point of line as pointObject.
            %
            %   Example: get start point of lineObject l1.
            %
            %       p = l1.startPoint();
            %
            p = obj.p1;
        end
        function p = endPoint(obj)
            %endPoint Get start point of line as pointObject
            %   p = endPoint() returns the end point of line as pointObject.
            %
            %   Example: get end point of lineObject l1.
            %
            %       p = l1.endPoint();
            %
            p = obj.p2;
        end
        function self = scale(self,scaleFactor)
            %scale Scale line dimensions and return line
            %   lineObject = scale(factor) returns a copy of current lineObject instance with its dimensions
            %   scaled by factor.
            %
            %   Example: scale dimensions of lineObject l1 by 25 and assign to l2.
            %
            %       l2 = l1.scale(25);
            %
            self.p1 = self.p1.scale(scaleFactor);
            self.p2 = self.p2.scale(scaleFactor);
            self.endParam = scaleFactor*self.endParam;
        end
        function value = isOn(self,p)
            %isOn Return true if point is on line, false otherwise
            %   isOn(p) returns true if pointObject p is on line and false otherwise.
            %
            %   Example: find if pointObject p1 is on lineObject l1.
            %
            %       isOn = l1.isOn(p1);
            %
            value = true;
            crossproduct = (p.y - self.p1.y) * (self.p2.x - self.p1.x) - (p.x - self.p1.x) * (self.p2.y - self.p1.y);
            lineLength = self.endParam;
            if abs(crossproduct) > eps*lineLength
                value = false;
            else
                dotproduct = (p.x - self.p1.x) * (self.p2.x - self.p1.x) + (p.y - self.p1.y)*(self.p2.y - self.p1.y);
                if dotproduct < 0
                    value = false;
                else
                    if dotproduct > lineLength*lineLength
                        value = false;
                    end
                end
            end
        end
    end
end