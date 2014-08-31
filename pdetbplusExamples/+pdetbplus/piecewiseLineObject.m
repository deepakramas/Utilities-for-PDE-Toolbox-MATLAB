classdef piecewiseLineObject < pdetbplus.boundaryObject
    %piecewiseLineObject Create piecewise line boundary from (x,y) table
    %   piecewiseLineObject is used for creating geometryObject instances.
    %
    %   piecewiseLineObject methods:
    %
    %       startPoint          - Get start point of piecewiseline as pointObject
    %       endPoint            - Get end point of piecewiseline as pointObject
    %       scale               - Scale piecewiseline dimensions and return piecewiseline
    %       translate           - Translate piecewiseline and return piecewiseline
    %       rotate              - Rotate piecewiseline with angle around origin and return piecewiseline
    %       getXY               - Return XY coordinates as [x,y], given distance along piecewiseline
    %       isOn                - Return true if point is on piecewiseline, false otherwise
    %
    %   piecewiseLineObject properties:
    %
    %       startParam          - Value of distance for start point of piecewiseline
    %       endParam            - Value of distance for end point of piecewiseline
    %       name                - Name of piecewiseline
    %       leftRegion          - Name of left region of piecewiseline going from start to end points
    %       rightRegion         - Name of right region of piecewiseline going from start to end points
    %
    %   Example: create piecewiseLineObject named 'A' out of segments (1,1)->(2,1.4)->(2.1,3)
    %
    %       l1 = piecewiseLineObject('A',[1,2,2.1],[1,1.4,3]);
    %
    %   See also boundaryObject arcObject lineObject parametricLineObject pointObject geometryObject.
    %
    % 
    properties (Hidden = true)
        x;
        y;
        params;
    end
    methods
        function self = piecewiseLineObject(name,x,y)
            %piecewiseLineObject Creates piecewise line from (x,y) table
            %   piecewiseLineObject(Name,X,Y) creates piecewiseLineObject named Name with x,y
            %   coordinates of line segments as X,Y.
            %
            %   Example: create piecewiseLineObject named 'A' out of segments (1,1)->(2,1.4)->(2.1,3)
            %
            %       l1 = piecewiseLineObject('A',[1,2,2.1],[1,1.4,3]);
            %
            self.name = name;
            self.x = x;
            self.y = y;
            self.startParam = 0;
            if length(x) ~= length(y)
                error 'unequal lengths for x and y';
            end
            if length(x) < 2
                error 'provide more than 1 point for x and y';
            end
            self.params = zeros(length(x),1);
            for k=2:length(x)
                dist = norm([x(k)-x(k-1),y(k)-y(k-1)]);
                self.params(k) = self.params(k-1) + dist;
            end
            [~,ix] = unique(self.params);
            self.params = self.params(ix);
            self.x = self.x(ix);
            self.y = self.y(ix);
            self.endParam = self.params(end);
        end
        function p = startPoint(self)
            import pdetbplus.*;
            [xpt,ypt] = self.getXY(self.startParam);
            p = pointObject(xpt,ypt);
        end
        function p = endPoint(self)
            import pdetbplus.*;
            [xpt,ypt] = self.getXY(self.endParam);
            p = pointObject(xpt,ypt);
        end
        function self = scale(self,scaleFactor)
            self.x = scaleFactor.*self.x;
            self.y = scaleFactor.*self.y;
            self.params = scaleFactor.*self.params;
            self.startParam = self.params(1);
            self.endParam = self.params(end);
        end
        function self = scaleX(self,scaleFactor)
            %scaleX Scale X dimensions of geometry            
            self.x = scaleFactor.*self.x;
            for k=2:length(self.x)
                dist = norm([self.x(k)-self.x(k-1),self.y(k)-self.y(k-1)]);
                self.params(k) = self.params(k-1) + dist;
            end
            self.startParam = self.params(1);
            self.endParam = self.params(end);
        end
        function self = scaleY(self,scaleFactor)
            %scaleY Scale Y dimensions of geometry
            self.y = scaleFactor.*self.y;
            for k=2:length(self.y)
                dist = norm([self.x(k)-self.x(k-1),self.y(k)-self.y(k-1)]);
                self.params(k) = self.params(k-1) + dist;
            end
            self.startParam = self.params(1);
            self.endParam = self.params(end);
        end        
        function self = translate(self,shiftXYPt)
            self.x = self.x + shiftXYPt.x;
            self.y = self.y + shiftXYPt.y;            
        end
        function self = rotate(self,origin,angle)
            import pdetbplus.*;
            for k=1:length(self.x)
                p = pointObject(self.x(k),self.y(k));
                p = p.rotate(origin,angle);
                self.x(k) = p.x;
                self.y(k) = p.y;
            end
        end
        function [xpt,ypt] = getXY(self,paramValue)
            paramValue = min(self.endParam,paramValue);
            xpt = zeros(length(paramValue),1);
            ypt = zeros(length(paramValue),1);
            for k=1:length(paramValue)
                indexUpper = find(self.params >= paramValue(k),1,'first');
                indexLower = find(self.params <= paramValue(k),1,'last');
                gridParams = [self.params(indexLower);self.params(indexUpper)];
                if (indexLower ~= indexUpper)
                    gridX = [self.x(indexLower);self.x(indexUpper)];
                    gridY = [self.y(indexLower);self.y(indexUpper)];
                    xpt(k) = interp1(gridParams,gridX,paramValue(k));
                    ypt(k) = interp1(gridParams,gridY,paramValue(k));
                else
                    xpt(k) = self.x(indexLower);
                    ypt(k) = self.y(indexLower);
                end
            end
        end        
        function [xmin,ymin,xmax,ymax] = getLimitsXY(self)
            %getLimitsXY Get maximum and minimum coordinates of boundary
            %   [xmin,ymin,xmax,ymax] = getLimitsXY returns x,y maximum and minimum coordinates
            %
            xmin = min(self.x);
            xmax = max(self.x);
            ymin = min(self.y);
            ymax = max(self.y);
        end
    end
end