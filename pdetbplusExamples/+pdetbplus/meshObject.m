classdef meshObject < pdetbplus.parseInputArgsObject
    %meshObject Create mesh
    %   meshObject creates mesh that is used by geometryObject. It calls initmesh from PDE Toolbox and
    %   therefore the mesh info in an meshObject instance can be used by other PDE Toolbox functions
    %   such as assempde, parabolic, hyperbolic.
    %
    %   meshObject methods:
    %       reset      - Reset mesh
    %       isEmpty    - Return true if mesh is empty, false otherwise
    %
    %   meshObject properties:
    %       p          - Point info for mesh as described in initmesh
    %       e          - Edge info for mesh as described in initmesh
    %       t          - Element info for mesh as described in initmesh
    %
    %   Example: create mesh from geometryObject g1 and set 'Hmax' param's value to 0.1,
    %       m1 = meshObject('geometry',g1.geometryFunction,'Hmax',0.1);
    %
    %   Example: create mesh from externally supplied mesh info : p,e,t
    %       m1 = meshObject('p',p,'e',e,'t',t);
    %
    %   See also geometryObject initmesh parabolic hyperbolic.
    %
    % 
    properties
        %p Point info for mesh as described in initmesh
        %   See also initmesh
        p;
        %e Edge info for mesh as described in initmesh
        %   See also initmesh
        e;
        %t Element info for mesh as described in initmesh
        %   See also initmesh
        t;
    end
    % fixed indices into e and t
    properties (Constant = true, Hidden = true)
        edgeLeftRegionIndex = 6;
        edgeRightRegionIndex = 7;
        edgeFirstPointIndex = 1;
        edgeSecondPointIndex = 2;
        triangleRegionIndex = 4;
        triangleFirstPointIndex = 1;
        triangleSecondPointIndex = 2;
        triangleThirdPointIndex = 3;
    end
    methods
        % instantiate either as
        % meshObject('p',p,'e',e,'t',t) or meshObject('geometry',h)
        function self = meshObject(varargin)
            %meshObject Creates mesh from geometry function handle or mesh nodes,edge,element info
            %   meshObject('geometry',h) creates mesh from geometry function handle h. 
            %   meshObject('p',p,'e',e,'t',t) creates mesh from mesh nodes,edge,element info of the format as
            %   described in the PDE Toolbox documentation.
            %   meshObject('numRefineMeshSteps',n) refines the mesh by calling PDE Toolbox function
            %   refinemesh n times and returns new mesh.
            %   Pass any name-value pairs that are valid for PDE Toolbox's initmesh function to the
            %   meshObject constructor.
            %
            %   Example: create mesh from geometryObject g1 and set 'Hmax' param's value to 0.1,
            %       m1 = meshObject('geometry',g1.geometryFunction,'Hmax',0.1);
            %
            %   Example: create mesh from externally supplied mesh info : p,e,t
            %       m1 = meshObject('p',p,'e',e,'t',t);
            %
            %   See also geometryObject initmesh.
            ps = self.parseInputArgs(varargin{:});           
            self.p = ps.Results.p;
            self.e = ps.Results.e;
            self.t = ps.Results.t;
            g = ps.Results.geometry;
            numRefineMeshSteps = ps.Results.numRefineMeshSteps;
            psu = fieldnames(ps.Unmatched);
            forwardedArgs = cell(0);
            for k=1:length(psu)
                forwardedArgs{end+1} = psu{k};
                forwardedArgs{end+1} = ps.Unmatched.(psu{k});
            end
            if ~isempty(g)
                self = self.init(g,numRefineMeshSteps,forwardedArgs{:});
            end
        end
        function self = reset(self)
            %reset Reset mesh
            self.p = [];
            self.e = [];
            self.t = [];
        end
        function unInitializedMesh = isEmpty(self)
            %isEmpty Return true if mesh is empty, false otherwise
            unInitializedMesh = isempty(self.p) || isempty(self.e) || isempty(self.t);
        end
    end
    methods (Hidden = true)
        function self = init(self, geometry, numRefineMesh,varargin)
            [self.p,self.e,self.t] = initmesh(geometry,varargin{:});
            for k=1:numRefineMesh
                [self.p,self.e,self.t]=refinemesh(geometry,self.p,self.e,self.t);
            end
            self.p = jigglemesh(self.p,self.e,self.t);
        end
    end
end