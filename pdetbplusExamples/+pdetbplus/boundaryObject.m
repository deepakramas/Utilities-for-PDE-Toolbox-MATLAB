classdef boundaryObject
    %boundaryObject Abstract class for objects that describe boundaries
    %   Boundary classes used in constructing geometryObject instances must subclass boundaryObject.
    %
    %   boundaryObject methods:
    %
    %       startPoint      - Get start point of boundary as pointObject
    %       endPoint        - Get end point of boundary as pointObject
    %       scale           - Scale boundary dimensions and return boundary
    %       translate       - Translate boundary and return boundary
    %       rotate          - Rotate boundary with angle around origin and return boundary
    %       getXY           - Return XY coordinates as [x,y], given param value along boundary
    %       isOn            - Return true if point is on boundary, false otherwise
    %
    %   boundaryObject properties:
    %
    %       startParam      - Value of parameter for start point of boundary
    %       endParam        - Value of parameter for end point of boundary
    %       name            - Name of boundary
    %       leftRegion      - Name of left region of boundary going from start to end points
    %       rightRegion     - Name of right region of boundary going from start to end points
    %
    %   See also geometryObject.
    %
    % 
    properties
        %startParam Value of parameter for start point of boundary
        startParam
        %endParam Value of parameter for end point of boundary
        endParam
        %name Name of boundary
        name
        %leftRegion Name of left region of boundary going from start to end points
        leftRegion
        %rightRegion Name of right region of boundary going from start to end points
        rightRegion
    end
    methods (Abstract)
        %startPoint Get start point of boundary as pointObject
        startPoint(obj)
        %endPoint Get end point of boundary as pointObject
        endPoint(obj)
        %scale Scale boundary dimensions and return boundary
        scale(self,scaleFactor)
        %translate Translate boundary and return boundary
        translate(self,shiftXYPt)
        %rotate Rotate boundary with angle around origin and return boundary
        rotate(self,origin,angle)
        %getXY Return XY coordinates as [x,y], given param value along boundary
        getXY(obj,paramValue)
    end
    methods
        function value = isOn(self,p)
            %isOn Return true if point is on boundary, false otherwise
            value = false;
            error 'isOn method undefined for boundary object\n';
        end
    end
end