classdef parametricLineObject < pdetbplus.boundaryObject
    %parametricLineObject Create parametrically described boundary
    %   parametricLineObject is used for creating geometryObject instances.
    %
    %   REQUIRES Symbolic Math Toolbox
    %
    %   parametricLineObject methods:
    %
    %       startPoint              - Get start point of boundary as pointObject
    %       endPoint                - Get end point of boundary as pointObject
    %       scale                   - Scale boundary dimensions and return boundary
    %       translate               - Translate boundary and return boundary
    %       rotate                  - Rotate boundary with angle around origin and return boundary
    %       getXY                   - Return XY coordinates as [x,y], given param value along boundary
    %       isOn                    - Return true if point is on boundary, false otherwise
    %
    %   parametricLineObject properties:
    %
    %       startParam              - Value of parameter for start point of boundary
    %       endParam                - Value of parameter for end point of boundary
    %       name                    - Name of boundary
    %       leftRegion              - Name of left region of boundary going from start to end points
    %       rightRegion             - Name of right region of boundary going from start to end points
    %
    %   Example: create segment named 'A' with parametric curve(2*cos(r),3*sin(r)) between r = 0 and
    %   r = pi
    %
    %       syms r;
    %       xSym = 2*cos(r);
    %       ySym = 3*sin(r);
    %       l1 = parametricLineObject('A',xSym,ySym,0,pi);
    %
    %   See also boundaryObject arcObject piecewiseLineObject parametricLineObject pointObject geometryObject.
    %
    %     
    properties (Hidden = true)
        xSym;
        ySym;
        paramSym;
        startLatch;
        endLatch;
    end
    methods
        function self = parametricLineObject(name,xExpression,yExpression,startParam,endParam,startLatch,endLatch)
            %parametricLineObject Create boundary segments from parametric curve 
            %   parametricLineObject(Name,xExpression,yExpression,startParam,endParam) creates
            %   a parametricLineObject named Name that has 
            %   coordinates of the currve as (x = xExpression(param), y = yExpression(param))
            %   startPoint = (xExpression(startParam), yExpression(startParam))
            %   endPoint = (xExpression(endParam), yExpression(endParam))
            %
            %   REQUIRES Symbolic Math Toolbox
            %
            %   Example: create segment named 'A' with parametric curve(2*cos(r),3*sin(r)) between r
            %   = 0 and r = pi
            %
            %       syms r;
            %       xSym = 2*cos(r);
            %       ySym = 3*sin(r); 
            %       l1 = parametricLineObject('A',xSym,ySym,0,pi);
            %
            self.name = name;
            self.xSym = sym(xExpression);
            self.ySym = sym(yExpression);
            xSymVar = symvar(self.xSym);
            ySymVar = symvar(self.ySym);
            if (xSymVar{1} ~= ySymVar{1})
                error 'xExpression and yExpression must have same parameter';
            end
            self.paramSym = xSymVar{1};
            self.startParam = startParam;
            self.endParam = endParam;
            if (nargin < 6)
                startLatch = [];
            end
            if (nargin < 7)
                endLatch = [];
            end
            self.startLatch = startLatch;
            self.endLatch = endLatch;
        end
        function p = startPoint(self)
            if isempty(self.startLatch)
                import pdetbplus.*;
                [x,y] = self.getXY(self.startParam);
                p = pointObject(x,y);
            else
                p = self.startLatch;
            end
        end
        function p = endPoint(self)
            if isempty(self.endLatch)
                import pdetbplus.*;
                [x,y] = self.getXY(self.endParam);
                p = pointObject(x,y);
            else
                p = self.endLatch;
            end
        end
        function self = scale(self,scaleFactor)
            self.xSym = self.xSym*scaleFactor;
            self.ySym = self.ySym*scaleFactor;
            if ~isempty(self.startLatch)
                self.startLatch = self.startLatch*scaleFactor;
            end
            if ~isempty(self.endLatch)
                self.endLatch = self.endLatch*scaleFactor;
            end
        end
        function self = translate(self,shiftXYPt)
            self.xSym = self.xSym + shiftXYPt.x;
            self.ySym = self.ySym + shiftXYPt.y;
            if ~isempty(self.startLatch)
                self.startLatch = self.startLatch + [shiftXYPt.x,shiftXYPt.y];
            end
            if ~isempty(self.endLatch)
                self.endLatch = self.endLatch + [shiftXYPt.x,shiftXYPt.y];
            end
        end
        function self = rotate(self,origin,angle)
            import pdetbplus.*;
            p = pointObject(self.xSym,self.ySym);
            p = p.rotate(origin,angle);
            self.xSym = p.x;
            self.ySym = p.y;
            if ~isempty(self.startLatch)
                self.startLatch = self.startLatch.rotate(origin,angle);
            end
            if ~isempty(self.endLatch)
                self.endLatch = self.endLatch.rotate(origin,angle);
            end
        end
        function [x,y] = getXY(self,paramValue)
            x = double(subs(self.xSym,self.paramSym,paramValue));
            y = double(subs(self.ySym,self.paramSym,paramValue));
            if ~isempty(self.startLatch) && ~isempty(find(paramValue == self.startParam,1))
                p = self.startPoint();
                x(paramValue == self.startParam) = p.x;
                y(paramValue == self.startParam) = p.y;
            elseif ~isempty(self.endLatch) && ~isempty(find(paramValue == self.endParam,1))
                p = self.endPoint();
                x(paramValue == self.endParam) = p.x;
                y(paramValue == self.endParam) = p.y;
            end
        end
    end
end