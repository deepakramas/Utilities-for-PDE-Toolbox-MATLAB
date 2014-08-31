classdef formulation < matlab.mixin.Copyable
    %formulation Superclass for equation forumulation aspect of PDE
    % 
    properties
        %geometry Geometry as geometryObject
        %   See also geometryObject
        geometry;
        %dimension Output dimension of PDE problem as int
        dimension = 1;
    end
    methods
        function self = formulation(g, dimension)
            %formulation Creates formulation
            %   formulation(g, dimension) creates formulation instance for geometryObject instance,
            %   g and output dimension of problem.
            %   See also geometryObject
            self.geometry = g;
            self.dimension = dimension;
        end
        function self = set.geometry(self,g)
            self.geometry = g;
        end
        function self = set.dimension(self,value)
            if value ~= floor(value) || value < 0
                error 'please enter positive integer value for dimension';
            end
            self.dimension = value;
        end
    end
end