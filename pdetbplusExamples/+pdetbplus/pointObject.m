classdef pointObject
    %pointObject  Represent point.
    %   pointObject represents a point (x,y) cartesian coordinates. It
    %   is used in the creation of geometryObject and other intermediate boundary
    %   related classes.
    %
    %   pointObject methods:
    %
    %       rotate                                  - Rotate point
    %       plus                                    - Add point
    %       minus                                   - Subtract point
    %       uplus                                   - Unary plus operator
    %       uminus                                  - Unary minus operator
    %       reflect                                 - Reflect point
    %       scale                                   - Scale point
    %       mtimes                                  - Multiply point
    %       mrdivide                                - Divide point
    %
    %   pointObject properties:
    %
    %       x                                       - X coordinate of point
    %       y                                       - Y coordinate of point
    %
    %   Example:  Create a point corresponding to (1,-2).
    %
    %      p1 = pointObject(1,-2);
    %  
    %   See also geometryObject.
    % 
    properties
        %x - X coordinate
        x;
        %y - Y coordinate
        y;
    end
    properties (Hidden = true)
        tags = cell(0);
    end
    methods
        function obj = pointObject(coords,coords2)
            %pointObject Create point
            %   pointObect(x,y) creates a point corresponding to (x,y)
            %
            %   Example: Create a point corresponding to (1,-2).
            %
            %       p1 = pointObject(1,-2);
            %
            if nargin < 2
                obj.x = coords(1);
                obj.y = coords(2);
            else
                obj.x = coords;
                obj.y = coords2;
            end
        end
        function obj = rotate(obj1,origin,angle)
            %rotate Rotate point anti-clockwise
            %   rotate(Origin,Angle) rotates point anti-clockwise around Origin with Angle
            %       Origin : Origin that point is rtoated about as pointObject
            %       Angle  : Angle of rotation in radians as double
            %
            %   Example: Rotate point, p1 around point (1,1),
            %   anti-clockwise by pi/4 radians and return rotated point, p2.
            %
            %       origin = pointObject(1,1);
            %       p2 = p1.rotate(origin,pi/4);
            %
            rotMat = [cos(angle) -sin(angle); sin(angle) cos(angle)];
            p = pdetbplus.pointObject(rotMat*[obj1.x - origin.x; obj1.y - origin.y]);
            obj = p + origin;
        end
        function obj = plus(obj1,obj2)
            %plus Add point
            %   + point2 adds point2 and returns point of type pointObject
            %       point2 can be of type pointObject OR an array of length 2 representing
            %       coordinates
            %
            %   Example: add point p2 of type pointObject to p1 of type pointObject and return p3 of
            %   type pointObject
            %
            %       p3 = p1 + p2;
            %
            %   Example: add [1,-1] to p1 of type pointObject and return p3 of type pointObject
            %
            %       p3 = p1 + [1,-1]; 
            %
            import pdetbplus.*;
            if isa(obj2,'pointObject')
                obj = pointObject([obj1.x+obj2.x;obj1.y+obj2.y]);
            else
                obj = pointObject([obj1.x+obj2(1);obj1.y+obj2(2)]);
            end
        end
        function obj = uplus(obj)
            %uplus Take unary positive
            %   + is a no-op and returns current point
            %
        end
        function obj = minus(obj1,obj2)
            %minus Subtract point
            %   - point2 subtracts point2 and returns point of type pointObject
            %       point2 can be of type pointObject OR an array of length 2 representing
            %       coordinates
            %
            %   Example: subtract point p2 of type pointObject to p1 of type pointObject and return p3 of
            %   type pointObject
            %
            %       p3 = p1 - p2;
            %
            %   Example: subtract [1,-1] to p1 of type pointObject and return p3 of type pointObject
            %
            %       p3 = p1 - [1,-1];
            %
            if isa(obj2,'pointObject')
                obj = pdetbplus.pointObject([obj1.x-obj2.x;obj1.y-obj2.y]);
            else
                obj = pdetbplus.pointObject([obj1.x-obj2(1);obj1.y-obj2(2)]);
            end
        end
        function obj = uminus(self)
            %uminus Take unary negative
            %   - negates the x,y values of current point and returns the negated value as a
            %   pointObject representing (-x,-y).
            %
            %   Example: return the negative of point p1 of type pointObject and assign it to p2
            %       p2 = -p1;
            obj = pointObject(-self.x,-self.y);
        end
        function obj = reflect(self,mirrorStartX,mirrorStartY,mirrorEndX,mirrorEndY)
            %reflect Reflect point
            %   reflect(x1,y1,x2,y2) reflects point across mirror represented by line with
            %   start point (x1,y1) and end point (x2,y2)
            %
            %   Example: reflect point, p1, across mirror line represented by start point (0,1) and end point
            %   (10,20) and assign result to point p2
            %
            %       p2 = p1.reflect(0,1,10,20);
            %
            xy = [self.x-mirrorStartX;self.y-mirrorStartY];
            mirror = [mirrorEndX-mirrorStartX;mirrorEndY-mirrorStartY];
            reflectXY = mirror*2*(mirror'*xy)/(mirror'*mirror) - xy;            
            obj = pdetbplus.pointObject(reflectXY(1)+mirrorStartX,reflectXY(2)+mirrorStartY);
        end
        function obj = scale(self,scaleFactor)
            %scale Scale by factor
            %   scale(factor) multiplies (x,y) coordinates of point by factor to return point with
            %   coordinates (x*factor,y*factor)
            %
            %   Example: scale coordinates of point p1 of type pointObject by factor 25 and return point p2 of type
            %   pointObject
            %
            %       p2 = p1.scale(25);
            %
            obj = pdetbplus.pointObject(scaleFactor*[self.x;self.y]);
        end
        function obj = mtimes(self,scaleFactor)
            %mtimes Multiply point by factor
            %   *factor scales (x,y) coordinates of point by factor to return point with coordinates
            %   (x*factor,y*factor). This function achieves an identical result as scale(factor).
            %
            %   Example: multiply coordinates of point p1 of type pointObject by 25 and return point p2 of type
            %   pointObject
            %
            %       p2 = p1 * 25;
            %
            %   See also scale
            obj = self.scale(scaleFactor);
        end
        function obj = mrdivide(self,scaleFactor)
            %mrdivide Divide point by factor
            %   /factor scales (x,y) coordinates of point by factor to return point with coordinates
            %   (x/factor,y*factor). This function achieves an identical result as scale(1/factor).
            %
            %   Example: divide coordinates of point p1 of type pointObject by 25 and return point p2 of type
            %   pointObject
            %
            %       p2 = p1 / 25;
            %
            %   See also scale
            obj = self.scale(1.0/scaleFactor);
        end        
    end
    methods (Hidden = true)
        function obj = addTag(tagName)
            obj.tags{end+1} = tagName;
        end
    end
end