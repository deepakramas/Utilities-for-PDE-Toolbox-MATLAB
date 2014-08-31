classdef geometryObject < pdetbplus.parseInputArgsObject
    %geometryObject Represent geometry and mesh for PDE Toolbox using boundary curves.
    %   geometryObject defines geometry by a set of bounding curves that are represented by
    %   subclasses of boundaryObject. A curve spits a geometry into left and right regions. In any
    %   geometry a curve should be defined exactly once. A curve should also not intersect other
    %   curves.
    %   A geometryObject instance can be meshed and passed to coeffsObject and
    %   boundaryConditionObject. Further a function handle can be passed as argument to PDE Toolbox
    %   function initmesh. More conveniently, you can pass the info in geometryObject's meshObject
    %   member, mesh, to assempde, parabolic and hyperbolic. geometryObject has modify, query, plot
    %   and post-processing methods.
    %
    %   geometryObject methods:
    %
    %       % Modify
    %       translate                           - Translate geometry
    %       scale                               - Scale geometry dimensions
    %       scaleX                              - Scale geometry along X dimension
    %       scaleY                              - Scale geometry along Y dimension
    %       rotate                              - Rotate geometry
    %       setInteriorExteriorRegions          - Assign interior and exterior region names
    %       plus                                - Add geometry
    %       initMesh                            - Initialize mesh
    %       % Plot
    %       plot                                - Plot geometry
    %       plotRegion                          - Plot Region
    %       % Query
    %       getRegionPoints                     - Get (x,y) coordinates of mesh points in region
    %       getRegionToId                       - Get numeric Id for region
    %       getRegionArea                       - Get area of region
    %       getTriangleAreasPerNode             - Get areas of triangles associated with mesh nodes
    %       getBoundaryEdges                    - Get boundary edges
    %       getBoundaryEdgesFromName            - Get boundary edges NOT IMPLEMENTED
    %       getRegionNodes                      - Get region nodes
    %       getBoundaryNodes                    - Get boundary nodes 
    %       getLimitsXY                         - Get maximum and minimum coordinates
    %       % Compute
    %       integrateOverRegion                 - Compute integral over region
    %       integrateOutsideBoundary            - Compute integral over boundary
    %       createXYFunctionFromNodalSolution   - Get (x,y) function of solution
    %       % Pre-defined shapes
    %       createCircle                        - Create circle
    %       createPolygon                       - Create polygon
    %       createSquare                        - Create square
    %       createGeometriesFromImage           - Create geometries from image
    %
    %   geometryObject properties:
    %
    %       name                                - Name of geometry of type char
    %       exteriorRegion                      - Name of exterior region
    %       mesh                                - Mesh of type meshObject
    %       geometryFunction                    - Function handle for geometry in PDE Toolbox as described in pdegeom
    %       regionToId                          - Map region names to Ids that are used by PDE Toolbox
    %
    %   Example: create triangle shaped geometryObject instance, g1, with name 'A' and cell array of
    %   subclasses of boundaryObject called b. The coordinates of g1 are (0,0), (1,1) and (1,-1).
    %   The inner region of g1 is 'tissue' and the outer region is 'fluid'. Plot and mesh g1.
    %
    %       p1 = pointObject(0,0); p2 = pointObject(1,1); p3 = pointObject(1,-1);
    %       b{1} = lineObject('l1',p1,p2); b{1}.leftRegion = 'fluid'; b{1}.rightRegion = 'tissue';
    %       b{2} = lineObject('l1',p2,p3); b{2}.leftRegion = 'fluid'; b{2}.rightRegion = 'tissue';
    %       b{3} = lineObject('l1',p3,p1); b{3}.leftRegion = 'fluid'; b{3}.rightRegion = 'tissue';
    %       g1 = geometryObject('A',b);
    %       g1.exteriorRegion = 'fluid';
    %       g1.plot();
    %       g1.initMesh('showMesh',true);
    %
    %   Example: use geometryObject instance g1 in assempde, adaptmesh
    %
    %       assempde(g1.geometryFunction,g1.mesh.p,g1.mesh.e,g1.mesh.t,...);
    %       adaptmesh(g1.geometryFunction,...);
    %
    %   See also boundaryObject lineObject arcObject piecewiseLineObject parametricLineObject
    %   pointObject coeffsObject boundaryConditionObject meshObject assempde adaptmesh parabolic hyperbolic.
    %
    %  
    properties (GetAccess = public, Hidden = false)
        %name Name of geometry of type char
        name;
        %exteriorRegion Name of exterior region. COMPULSORY to set prior to meshing or plotting.
        exteriorRegion;
        %mesh Mesh of type meshObject
        %   See also meshObject
        mesh;
    end
    properties (SetAccess = private)
        %boundary List of boundaries that define the mesh
        boundary = cell(0);
        %regions Regions or "sub-domains" of the geometry
        regions = cell(0);
    end
    properties (SetAccess = private, Hidden = true)
        dlCache;
    end
    properties (Dependent = true,SetAccess = private)
        %geometryFunction Function handle for geometry in PDE Toolbox as described in pdegeom
        %   geometryFunction is an argument to assempde(), adaptmesh(), parabolic() and hyperbolic()
        %   in PDE Toolbox. This handle can also be passed to initmesh() from PDE Toolbox.
        %
        %   See also pdegeom initmesh assempde adaptmesh parabolic hyperbolic. 
        geometryFunction;
        %regionToId Map region names to Ids that are used by PDE Toolbox
        regionToId;
    end
    methods
        function self = geometryObject(name,c)
            %geometryObject Create geometry
            %   geometryObject(Name,BoundaryCellArray) creates geometry with name equal to Name and
            %   a cell array of subclasses of boundaryObject.
            %
            %   Example: create geometryObject instance, g1, corresponding to name 'A' and cell array of
            %   subclasses of boundaryObject called b.
            %
            %       g1 = geometryObject('A',b);
            %
            %   See also boundaryObject.
            import pdetbplus.meshObject;
            for k=1:length(c)
                self.boundary{k} = c{k};
                self.regions{end+1} = c{k}.leftRegion;
                self.regions{end+1} = c{k}.rightRegion;
            end
            self.regions = unique(self.regions);
            self.name = name;
            self.mesh = meshObject();
        end
        function gf = get.geometryFunction(self)
            self = self.createBoundaryFunctionCache(); % set cache to avoid repeated calls
            gf = @self.geometryFunction_impl;
        end
        function exr = get.exteriorRegion(self)
            if isempty(self.exteriorRegion)
                error 'exterior region must be defined prior to plotting or meshing';
            end
            exr = self.exteriorRegion;
        end
        function regionToId = get.regionToId(self)
            regionToId = containers.Map;
            ix = 0;
            regionToId(self.exteriorRegion) = ix;
            for k=1:length(self.regions)
                if ~strcmp(self.regions{k},self.exteriorRegion)
                    ix = ix + 1;
                    regionToId(self.regions{k}) = ix;
                end
            end
        end
        function self = translate(self,t)
            %translate Translate geometry
            %   geometryObject = translate(p) translates geometry by direction specified by
            %   pointObject p and returns the a new translated geometryObject instance.
            %
            %   Example: translate geometry g1 by (1,-3) and assign to new geometry g2.
            %       p = pointObject(1,-3);
            %       g2 = g1.translate(p);
            %
            %   See also pointObject.
            for k=1:length(self.boundary)
                self.boundary{k} = self.boundary{k}.translate(t);
            end
            if ~self.mesh.isEmpty()
                self = self.initMesh();
            end
        end
        function self = scale(self,scaleFactor)
            %scale Scale geometry dimensions
            %   geometryObject = scale(factor) scales geometry dimensions by factor and returns a
            %   new scaled geometryObject instance.
            %
            %   Example: scale geometry g1 by 25 and assign to new geometry g2.
            %       g2 = g1.scale(25);
            %
            for k=1:length(self.boundary)
                self.boundary{k} = self.boundary{k}.scale(scaleFactor);
            end
        end
        function self = scaleX(self,scaleFactor)
            %scaleX Scale geometry along X dimension
            %   geometryObject = scaleX(factor) scales geometry along X dimension by factor and returns a
            %   new scaled geometryObject instance.
            %
            %   Example: scale X dimension of geometry g1 by 25 and assign to new geometry g2.
            %       g2 = g1.scaleX(25);
            %
            %   See also scaleY.
            for k=1:length(self.boundary)
                b = self.boundary{k};
                if ismethod(b,'scaleX')
                    self.boundary{k} = self.boundary{k}.scaleX(scaleFactor);
                end
            end
        end
        function self = scaleY(self,scaleFactor)
            %scaleY Scale geometry along Y dimension
            %   geometryObject = scaleY(factor) scales geometry along Y dimension by factor and returns a
            %   new scaled geometryObject instance.
            %
            %   Example: scale Y dimension of geometry g1 by 25 and assign to new geometry g2.
            %       g2 = g1.scaleY(25);
            %
            %   See also scaleX.
            for k=1:length(self.boundary)
                b = self.boundary{k};
                if ismethod(b,'scaleY')
                    self.boundary{k} = self.boundary{k}.scaleY(scaleFactor);
                end
            end
        end        
        function self = rotate(self,origin,angle)
            %rotate Rotate geometry
            %   geometryObject = rotate(Origin,Angle) rotates geometry by Angle radians
            %   anti-clockwise around Origin that is of type pointObject and returns a new
            %   geometryObject instance.
            %
            %   Example: rotate geometry g1 by pi/4 around (2,3) and assign to new geometry g2.
            %       origin = pointObject(2,3);
            %       g2 = g1.rotate(origin,pi/4);
            %
            %   See also pointObject.
            for k=1:length(self.boundary)
                self.boundary{k} = self.boundary{k}.rotate(origin,angle);
            end
            if ~self.mesh.isEmpty()
                self = self.initMesh();
            end
        end
        function self = setInteriorExteriorRegions(self,interiorRegionName,exteriorRegionName,isClockwise)
            %setInteriorExteriorRegions Assign interior and exterior region names
            %   geometryObject =
            %   setInteriorExteriorRegions(interiorRegionName,exteriorRegionName,isClockwise) sets
            %   the names of interior and exterior regions of a geometry. The orientation i.e.
            %   clockwise or anti-clockwise of the boundaries is also a required input. Do NOT use
            %   this method unless the geometry is closed and orientation is known.
            %
            %   Example: See createGeometriesFromImage for an example.
            %
            %   See also createGeometriesFromImage.
            if isClockwise
                leftRegionName = exteriorRegionName;
                rightRegionName = interiorRegionName;
            else
                rightRegionName = exteriorRegionName;
                leftRegionName = interiorRegionName;
            end
            for k=1:length(self.boundary)
                self.boundary{k}.leftRegion = leftRegionName;
                self.boundary{k}.rightRegion = rightRegionName;
            end
            self.regions = cell(2,1);
            self.regions{1} = leftRegionName;
            self.regions{2} = rightRegionName;
            self.exteriorRegion = exteriorRegionName;
        end
        function self = plus(self,obj2)
            %plus Add geometry
            %   + g2 adds geometryObject, g2 to current geometryObject instance and returns new
            %   geometryObject instance. The boundaries of g2 should NOT intersect the boundaries of
            %   current geometryObject instance.
            %
            %   Example: add geometryObject, g2 to geometryObject, g1 and assign it to
            %   geometryObject, g3.
            %       g3 = g1 + g2;
            self = self.add(obj2);
        end
        function self = initMesh(self,varargin)
            %initMesh Initialize mesh
            %   geometryObject = initMesh(Param1,Value1,Param2,Value2,...) initializes mesh
            %   corresponding to this geometry where valid values for PARAM are:
            %       'showMesh' : Flag for displaying mesh as bool. Default is false.
            %       Any Param-Value pair that is valid for PDE Toolbox's initmesh() can
            %       be passed.
            %   A new instance of geometryObject that has an initialized mesh is returned.
            %
            %   Example: initialize mesh for geometry g1 displaying mesh and setting 'Hmax' to be
            %   0.1 and assign new geometryObject instance to g2.
            %       g2 = g1.initMesh('showMesh',true,'Hmax',0.1);
            %
            %   See also initmesh plot meshObject.
            import pdetbplus.meshObject;
            p = self.parseInputArgs(varargin{:});
            pu = fieldnames(p.Unmatched);
            forwardedArgs = cell(0);
            for k=1:length(pu)
                forwardedArgs{end+1} = pu{k};
                forwardedArgs{end+1} = p.Unmatched.(pu{k});
            end            
            h = self.geometryFunction;
            self.mesh = meshObject('geometry',h,'numRefineMeshSteps',p.Results.numRefineMeshSteps,forwardedArgs{:});
            if nargin > 1 && p.Results.showMesh
                self.plot('showMesh',p.Results.showMesh);
            end
        end
        function [] = plot(self,varargin)
            %plot Plot geometry
            %   plot(Param,Value) where valid values for PARAM are:
            %       'showMesh' : Flag for displaying mesh in plot as bool. Default is false. Mesh is
            %       should be initialized if 'showMesh' is set to true.
            %
            %   Example: plot geometry g1 and display mesh.
            %       g1.plot('showMesh',true);
            %
            %   See also plotRegion initMesh.
            p = self.parseInputArgs(varargin{:});
            if ~p.Results.showMesh                
                h = self.geometryFunction;                
                pdegplot(h);
            else
                if isempty(self.mesh) || isempty(self.mesh.p)
                    error 'call method initMesh() prior to calling plot';
                else
                   pdemesh(self.mesh.p,self.mesh.e,self.mesh.t); 
                end
            end
        end
        function [] = plotRegion(self,regionName)
            %plotRegion Plot region
            %   plotRegion(regionName) plots region named regionName. The current geometryObject
            %   instance MUST be meshed prior to calling plotRegion.
            %
            %   Example: plot region 'A' in geometry g1.
            %       g1.plotRegion('A');
            %
            %   See also plot initMesh.
            self.plot();
            hold on;
            if isempty(self.mesh)
                error 'mesh must be initialized before plotting individual region';
            end
            tRegion = self.mesh.t(1:3,self.mesh.t(self.mesh.triangleRegionIndex,:) == self.regionToId(regionName));
            if ~isempty(tRegion)
                tr = TriRep(tRegion',self.mesh.p');
                triplot(tr);
                title(regionName);
            end
        end
        function [xpts,ypts] = getRegionPoints(self,regionName)
            %getRegionPoints Get (x,y) coordinates of mesh points in region.
            %   [xpts,ypts] = getRegionPoints(regionName) returns (x,y) coordinates of all mesh
            %   nodes in region named regionName in geometry.
            %
            %   Example: get (x,y) coordinates for mesh nodes in region 'A' in geometry g1.
            %       [xpts,ypts] = g1.getRegionPoints('A');
            %
            %   See also getBoundaryEdges getBoundaryEdgesFromName getRegionNodes getBoundayNodes.
            pts = cell(0);
            for k=1:length(self.boundary)
                b = self.boundary{k};
                if strcmp(b.leftRegion,regionName) || strcmp(b.rightRegion,regionName)
                    pts{end+1} = b.startPoint;
                    pts{end+1} = b.endPoint;
                end
            end
            xpts = zeros(length(pts),1);
            ypts = zeros(length(xpts),1);
            for k=1:length(pts)
                xpts(k) = pts{k}.x;
                ypts(k) = pts{k}.y;
            end
        end      
        function regionId = getRegionToId(self,regionName)
            %getRegionToId Get numeric Id for region
            %   regionId = getRegionToId(regionName) returns numeric Id, regionId, for region named
            %   regionName. The id can be used in boundary condition and coefficient functions of
            %   PDE Toolbox.
            %
            %   Example: get numeric Id for region called 'A' in geometry g1.
            %       id = g1.getRegionToId('A');
            %
            regionToIds = self.regionToId;
            if isempty(regionToIds)
                error 'region names have not been mapped to region Ids';
            end
            regionId = regionToIds(regionName);
        end
        function area = getRegionArea(self,regionName)
            %getRegionArea Get area of region
            %   area = getRegionArea(regionName) returns area of region named regionName.
            %
            %   Example: get area of region called 'A' in geometry g1.
            %       areaA = g1.getRegionArea('A');
            %
            %   See also getTriangleAreasPerNode.
            regionToIds = self.regionToId;
            [ar,~] = pdetrg(self.mesh.p,self.mesh.t(:,self.mesh.t(4,:) == regionToIds(regionName)));
            area = sum(ar);
        end
        function areasPerNode = getTriangleAreasPerNode(self,regionName)
            %getTriangleAreasPerNode Get areas of triangles associated with mesh nodes
            %   areasPerNode = getTriangleAreasPerNode(regionName) returns cell array of areas
            %   attached to each mesh node. The length of areasPerNode is #mesh nodes.
            %
            %   Example: get cell array of areas per node for region 'A' in geometry g1.
            %       areasPerNode = g1.getTriangleAreasPerNode('A');
            %
            %   See also getRegionArea.   
            regionToIds = self.regionToId;
            t = self.mesh.t(:,self.mesh.t(4,:) == regionToIds(regionName));
            [ar,~] = pdetrg(self.mesh.p,t);
            areasPerNode = cell(size(self.mesh.p,2),1);
            for k=1:length(ar)
                areasPerNode{t(1,k)} = [areasPerNode{t(1,k)} ar(k)];
                areasPerNode{t(2,k)} = [areasPerNode{t(2,k)} ar(k)];
                areasPerNode{t(3,k)} = [areasPerNode{t(3,k)} ar(k)];
            end
        end
        function [boundaryEdges,boundaryEdgeElements] = getBoundaryEdges(self,regionAName,regionBName)
            %getBoundaryEdges Get boundary edges
            %   boundaryEdges = getBoundaryEdges(regionA,regionB) returns mesh edges for boundary
            %   between regions named regionA and regionB.
            %
            %   Example: get boundary mesh edges between regions 'A' and 'B' for geometry g1.
            %       boundaryEdges = g1.getBoundaryEdges('A','B');
            %
            %   See also getBoundaryEdgesFromName getBoundaryEdgesFromName getRegionNodes getRegionPoints.
            regionAId = self.regionToId(regionAName);
            regionBId = self.regionToId(regionBName);
            boundaryEdges = union(intersect(find(self.mesh.e(self.mesh.edgeLeftRegionIndex,:) == regionAId), find(self.mesh.e(self.mesh.edgeRightRegionIndex,:) == regionBId)),...
                intersect(find(self.mesh.e(self.mesh.edgeLeftRegionIndex,:) == regionBId), find(self.mesh.e(self.mesh.edgeRightRegionIndex,:) == regionAId)));
            if nargout > 1
                pMatrix = sparse(size(self.mesh.p,2),size(self.mesh.p,2));
                for k=1:size(self.mesh.t,2)
                    p1 = self.mesh.t(1,k);
                    p2 = self.mesh.t(2,k);
                    p3 = self.mesh.t(3,k);
                    pMatrix([p1 p2 p3],[p1 p2 p3]) = k*ones(3,3);
                end
                boundaryEdgeElements = zeros(size(boundaryEdges));
                count = 0;
                for edge=boundaryEdges
                    count = count + 1;
                    p1 = self.mesh.e(1,edge);
                    p2 = self.mesh.e(2,edge);
                    boundaryEdgeElements(count) = pMatrix(p1,p2);
                end
                warning off;
                tr = TriRep(self.mesh.t(1:3,boundaryEdgeElements)',self.mesh.p');
                warning on;
                triplot(tr);
            end
        end
        function boundaryEdges = getBoundaryEdgesFromName(self,boundaryName)
            %getBoundaryEdgesFromName Get boundary edges NOT IMPLEMENTED
            %
            %   See also getBoundaryEdgesFromName getBoundaryEdges getRegionNodes getRegionPoints.
            error 'this function has not been implemented';
        end        
        function regionNodes = getRegionNodes(self,varargin)
            %getRegionNodes Get region nodes
            %   regionNodes = getRegionNodes('name', regionName) returns all indices for mesh nodes
            %   on region called regionName.
            %
            %   Example: get region mesh nodes of region called 'A' in geometry g1.
            %       regionNodes = g1.getRegionNodes('name','A');
            %
            %   See also getBoundaryNodes getBoundaryEdges getBoundaryEdgesFromName getRegionPoints..            
            p = self.parseInputArgs(varargin{:});
            regionName = p.Results.name;
            regionNodes = [];
            hold on;
            if isempty(self.mesh)
                error 'mesh must be initialized before plotting individual region';
            end
            tRegion = self.mesh.t(1:3,self.mesh.t(self.mesh.triangleRegionIndex,:) == self.regionToId(regionName));
            if ~isempty(tRegion)
                regionNodes = unique(sort([tRegion(1,:) tRegion(2,:) tRegion(3,:)]));
            end
        end
        function boundaryNodes = getBoundaryNodes(self,varargin)
            %getBoundaryNodes Get boundary nodes
            %   boundaryNodes = getBoundaryNodes(Param1, Value1, ...) returns all indices for mesh nodes
            %   on boundary. Valid values for PARAM are:
            %       'name'        : Name of boundary as char. Optional.
            %       'region1'     : Name of first region defining boundary as char. Optional.
            %       'region2'     : Name of second region defining boundary as char. Optional.
            %   Either 'name' OR both regions are input.
            %
            %   Example: get boundary nodes of boundary between regions 'A' and 'B' for geometry g1.
            %   Boundary is called 'b1'.
            %       boundaryNodes = g1.getBoundaryNodes('region1','A','region2','B');
            %       % OR use boundary name syntax
            %       boundaryNodes = g1.getBoundaryNodes('name','b1');
            %
            %   See also getRegionNodes getBoundaryEdges getBoundaryEdgesFromName getRegionPoints.
            p = self.parseInputArgs(varargin{:});
            region1Name = p.Results.region1;
            region2Name = p.Results.region2;
            boundaryName = p.Results.name;
            boundaryNodes = [];
            if ~isempty(region1Name) && ~isempty(region2Name)
                boundaryEdges = self.getBoundaryEdges(region1Name,region2Name);
                boundaryNodes = unique([self.mesh.e(1,boundaryEdges) self.mesh.e(2,boundaryEdges)]);
            elseif ~isempty(boundaryName)
                for bcix=1:length(self.boundary)
                    if strcmp(self.boundary{bcix}.name,boundaryName) && ismethod(self.boundary{bcix},'isOn')
                        for pix = 1:size(self.mesh.p,2)
                            if self.boundary{bcix}.isOn(pdetbplus.pointObject(self.mesh.p(1,pix),self.mesh.p(2,pix)))
                                boundaryNodes = [boundaryNodes pix];
                            end
                        end
                    end
                end
            end
        end
        function [xmin,ymin,xmax,ymax] = getLimitsXY(self)
            %getLimitsXY Get maximum and minimum coordinates
            %   [xmin,ymin,xmax,ymax] = getLimitsXY returns the coordinates limits of geometry
            %       xmin : Minimum X dimension of geometry as double
            %       xmax : Maximum X dimension of geometry as double
            %       ymin : Minimum Y dimension of geometry as double
            %       ymax : Maximum Y dimension of geometry as double
            %
            %   Example: get maximum and minimum coordinates amongst all points in geometry, g1.
            %       [xmin,ymin,xmax,ymax] = g1.getLimitsXY();
            %
            xmin = realmax;
            ymin = xmin;
            xmax = realmin;
            ymax = xmax;
            for k=1:length(self.boundary)
                b = self.boundary{k};
                if ismethod(b,'getLimitsXY')
                    [xmintmp,ymintmp,xmaxtmp,ymaxtmp] = b.getLimitsXY();
                    xmin = min(xmintmp,xmin);
                    ymin = min(ymintmp,ymin);
                    xmax = max(xmaxtmp,xmax);
                    ymax = max(ymaxtmp,ymax);
                else
                    warning 'getLimitsXY not implemented for a boundary type of region; returning null';
                    xmin = [];
                    ymin = xmin;
                    xmax = xmin;
                    ymax = xmin;
                    break;
                end
            end
        end
        function value = integrateOverRegion(self,varargin)
            %integrateOverRegion Compute integral over region
            %   value = integrateOverRegion(param1,Value1,...) returns the value of an
            %   integration computed over region. Valid values for PARAM are:
            %       'region'            : Name of region over which integration is performed as
            %       char.
            %       'integrand'         : Integrand function as function handle. The inputs to the
            %       integrand function are (xmid,ymid,additionalInputs) where (xmid,ymid) is the
            %       midpoint of a triangle in the region. xmid,ymid are automatically input by
            %       integrateOverRegion. additionalInputs are additional input as defined via
            %       'integrandInput'.
            %       'integrandInput'    : Additional input arguments for integrand as cell array.
            %       Each element of cell array is an array defined over the mesh triangles and is of
            %       size length(#mesh triangles) X 1.
            %
            %   Example: Integrate function integrandSegs over region 'fooMaterial' for geometry
            %   g1.
            %
            %       function val = integrandSegs(xmid,ymid,inputs)
            %           umid = inputs{1};
            %           val = umid*sqrt(xmid^2+ymid^2); % some function
            %       end
            %       % Get solution over mesh triangles; u is solution over nodes.
            %       utriangles = pdeintrp(g1.mesh.p,g1.mesh.t,u);
            %       val1 = g1.integrateOverRegion('integrand',@integrandSegs,'region','fooMaterial','integrandInput',{utriangles});
            %
            %   See also integrateOutsideBoundary mesh.
            ps = self.parseInputArgs(varargin{:});
            regionName = ps.Results.region;
            f = ps.Results.integrand;
            additionalInputs = ps.Results.integrandInput;
            regionToIds = self.regionToId;
            regionId = regionToIds(regionName);
            regionIndices = find(self.mesh.t(self.mesh.triangleRegionIndex,:) == regionId);
            ar = pdetrg(self.mesh.p,self.mesh.t(:,regionIndices));
            value = 0;
            for k = 1: length(regionIndices)
                t = self.mesh.t(1:3,regionIndices(k));
                xymid = sum(self.mesh.p(1:2,t),2)/3;
                additionalInputsAtTriangle = cell(length(additionalInputs),1);
                for kx=1:length(additionalInputs)
                    tmpInput = additionalInputs{kx}(regionIndices(k));
                    additionalInputsAtTriangle{kx} = tmpInput;
                end
                fvalue = f(xymid(1),xymid(2),additionalInputsAtTriangle);
                value = value + fvalue*ar(k);
            end
        end
        function value = integrateOutsideBoundary(self,varargin)
            %integrateOutsideBoundary Compute integral over boundary
            %   value = integrateOutsideBoundary(param1,Value1,...) returns the value of an
            %   integration computed just outside the specified boundary. Valid values for PARAM are:
            %       'integrand'         : Integrand function as function handle. The inputs to the
            %       integrand function are (xmid,ymid,nxmid,nymid,additionalInputs) where
            %       (xmid,ymid) is the midpoint of an edge on the boundary, (nxmid,nymid) is the
            %       normal on the edge. xmid,ymid,nxmid,nymid are automatically input by
            %       integrateOutsideBoundary. additionalInputs are additional input as defined via
            %       'integrandInput'.
            %       'integrandInput'    : Additional input arguments for integrand as cell array.
            %       Each element of cell array is an array defined over the mesh triangles and is of
            %       size length(#mesh triangles) X 1.
            %       'doNotSum'          : Flag for not summing the values of the integrals on the
            %       boundary edges as bool. Default is false. If set to false, the integral over all
            %       boundary edges is returned. If set to true, an #mesh edges X 1 array is returned
            %       with value only set at edges on the boundary.
            %       'name'              : Name of boundary as char. Optional. Normals for edges on the
            %       boundary point from the RHS side to the LHS side of the edges. Either 'name' OR
            %       both regions are input.
            %       'innerRegion'       : Name of inner region defining boundary as char. Optional.
            %       Boundary normal points from inner region to outer region. Either 'name' OR both
            %       regions are input.
            %       'outerRegion'       : Name of outer region defining boundary as char. Optional.
            %       Boundary normal points from inner region to outer region. Either 'name' OR both
            %       regions are input.
            %       'showBoundary'      : Plot the boundary mesh triangles corresponding to the
            %       integration as bool. Default is true. Plotted normals are scaled to edge
            %       lengths.
            %       'showBoundaryClear' : Clear current axes prior to showing boundary. Default is
            %       false.
            %   If 'integrandInputs' is specified, 'outerRegion' either explicitly specified or
            %   implied by use of 'name' MUST be meshed.
            %
            %   Example: Integrate 'integrandNeumann' function on boundary defined by inner region =
            %   'air' and outer region = 'fooMaterial' in geometry g1.
            %
            %       function val = integrandNeumann(xmid,ymid,nxmid,nymid,inputs)
            %`          % xmid,ymid,nxmid,nymid are automatically input by integrateOutsideBoundary.
            %           ux = inputs{1};
            %           uy = inputs{2};
            %           val = (ux*nxmid + uy*nymid)^2;
            %       end
            %       % u is the solution, (cdudx,cdudy) are the C.grad components defined over all mesh triangles 
            %       [cdudx,cdudy] = pdecgrad(g1.mesh.p,g1.mesh.t,coeffs.cFunction,u);
            %       val1 = g1.integrateOutsideBoundary('innerRegion','air','outerRegion','fooMaterial','integrand',@integrandNeumann,...
            %           'showBoundary',true,'integrandInput',{cdudx,cdudy});
            %
            %   See also integrateOverRegion.
            ps = self.parseInputArgs(varargin{:});
            outerRegion = ps.Results.outerRegion;
            innerRegion = ps.Results.innerRegion;
            boundaryName = ps.Results.name;
            showBoundary = ps.Results.showBoundary;
            showBoundaryClear = ps.Results.showBoundaryClear; % esoteric but needed
            doNotSum = ps.Results.doNotSum;
            f = ps.Results.integrand;
            integrandInputs = ps.Results.integrandInput;
            value = 0;
            regionToIds = self.regionToId;
            if isempty(boundaryName)
                % find out boundary corresponding to inner and outer regions
                outerRegionId = regionToIds(outerRegion);
                innerRegionId = regionToIds(innerRegion);
                regionIndices = find(self.mesh.t(self.mesh.triangleRegionIndex,:) == outerRegionId);
                regionTriangles = self.mesh.t(1:3,regionIndices);
                rt{1} = regionTriangles(1:2,:);
                rt{2} = regionTriangles(2:3,:);
                rt{3} = regionTriangles(3:-2:1,:);
                rt{4} = regionTriangles(2:-1:1,:);
                rt{5} = regionTriangles(3:-1:2,:);
                rt{6} = regionTriangles(1:2:3,:);
                additionalInputs = cell(length(integrandInputs),1);
                for kx=1:length(integrandInputs)
                    additionalInputs{kx} = integrandInputs{kx}(regionIndices);
                end
                boundaryEdges = union(intersect(find(self.mesh.e(self.mesh.edgeLeftRegionIndex,:) == innerRegionId), find(self.mesh.e(self.mesh.edgeRightRegionIndex,:) == outerRegionId)),...
                    intersect(find(self.mesh.e(self.mesh.edgeLeftRegionIndex,:) == outerRegionId), find(self.mesh.e(self.mesh.edgeRightRegionIndex,:) == innerRegionId)));
            else
                boundaryEdges = self.getBoundaryEdgesFromName(boundaryName);
            end
            borderTrianglesIndex = [];
            Nx1 = []; Nx2 = []; Ny1 = []; Ny2 = [];
            for ix=1:length(boundaryEdges)
                p1 = self.mesh.e(self.mesh.edgeFirstPointIndex,boundaryEdges(ix));
                p2 = self.mesh.e(self.mesh.edgeSecondPointIndex,boundaryEdges(ix));
                xy1 = self.mesh.p(1:2,p1);
                xy2 = self.mesh.p(1:2,p2);
                dp = xy2 - xy1;
                lengthEdge = norm(dp);
                nxe = -dp(2)/lengthEdge;
                nye = dp(1)/lengthEdge;
                if (self.mesh.e(self.mesh.edgeLeftRegionIndex,boundaryEdges(ix)) == innerRegionId) && (self.mesh.e(self.mesh.edgeRightRegionIndex,boundaryEdges(ix)) == outerRegionId)
                    % nxe, nye points from outer region to inner region. we need to flip the direction
                    nxe = -nxe;
                    nye = -nye;
                end
                xymid = (xy1 + xy2)/2;
                % find corresponding mesh element i.e. triangle
                triangleIndex = [];
                for kx=1:length(rt)
                    triangleIndex = find((rt{kx}(1,:)==p1) & (rt{kx}(2,:)==p2));
                    if triangleIndex
                        borderTrianglesIndex = [borderTrianglesIndex triangleIndex];
                        Nx1 = [Nx1 xymid(1)];
                        Nx2 = [Nx2 (nxe*lengthEdge)];
                        Ny1 = [Ny1 xymid(2)];
                        Ny2 = [Ny2 (nye*lengthEdge)];
                        break;
                    end
                end
                if triangleIndex
                    additionalInputsAtTriangle = cell(length(additionalInputs),1);
                    for kx=1:length(additionalInputs)
                        additionalInputsAtTriangle{kx} = additionalInputs{kx}(triangleIndex);
                    end
                    fvalue = f(xymid(1),xymid(2),nxe,nye,additionalInputsAtTriangle);
                    if ~doNotSum
                        value = value + fvalue*lengthEdge;
                    else
                        value(boundaryEdges(ix),1) = fvalue*lengthEdge;
                    end
                end
            end
            if showBoundary && ~isempty(regionTriangles)
                warning off; % creates warning that is not relevant here
                tr = TriRep(regionTriangles(:,borderTrianglesIndex)',self.mesh.p');
                warning on;
                if showBoundaryClear
                    cla;
                end
                triplot(tr);
                hold on;
                quiverHandle = quiver(Nx1,Ny1,Nx2,Ny2);
                set(quiverHandle,'Color','red');
            end      
        end
        function xyFunctionInstance = createXYFunctionFromNodalSolution(self,sol)
            %createXYFunctionFromNodalSolution Get (x,y) function of solution
            %   Fxy = createXYFunctionFromNodalSolution(Solution) creates a function handle from the
            %   solution of a PDE that has been defined for the geometry
            %       Solution : Solution for a PDE that has been defined on the geometry as an N x M
            %       array where N = # nodes, M = #time points. Solution must be for the mesh of
            %       this geometryObject.
            %       Fxy      : Function handle that returns a max(length(x),length(y)) X
            %       size(solution,2) array for any coordinate input (x,y) within this
            %       geometryObject. For any point outside geometryBoject, the value is NaN
            %
            %   Example: Return function handle f for solution sol and geometry g1
            %       f = g1.createXYFunctionFromNodalSolution(sol);
            %       % Call f on point (1,-1) within g1
            %       fval = f(1,-1);
            %       % Call f on points (1,1) and (1,-1) within g1
            %       fvals = f([1,1],[1,-1]);
            %       % Another way to call
            %       fvals = f(1,[1,-1]);
            %
            %   See also assempde
            xyFunc = pdetbplus.xyFunctionClass(sol,self.mesh);
            xyFunctionInstance = @xyFunc.xyFunction;
        end
    end
    methods(Hidden = true)
        function self = add(self,obj2)      
            for k=1:length(obj2.boundary)
                self.boundary{end+1} = obj2.boundary{k};
            end
            for k=1:length(obj2.regions)
                self.regions{end+1} = obj2.regions{k};
            end
            self.regions = unique(self.regions);
        end
    end
    %% Internal helper methods
    methods(Access = private, Hidden = true)
        function [x,y] = geometryFunction_impl(self,bs,s)
            %method for returning geometry function. pass handle to this to PDE tb functions. uses
            %boundary rep to create geometry
            
            nbs = length(self.boundary);
            if nargin == 1
                x = nbs;
                return;
            end
            dl = zeros(4,length(self.boundary));
            if ~isempty(self.dlCache)
                dl = self.dlCache;
            else
                regionToId = self.regionToId;
                for k=1:length(self.boundary)
                    b = self.boundary{k};
                    leftRegionId = regionToId(b.leftRegion);
                    rightRegionId = regionToId(b.rightRegion);
                    dl(:,k) = [b.startParam b.endParam leftRegionId rightRegionId]';
                end
            end
            
            if nargin==2
                x=dl(:,bs);
                return
            end
            transposed = false;
            if size(s,2) == 1
                s = s.';
                bs = bs.';
                transposed = true;
            end
            
            [m,n]=size(bs);
            if m==1 && n==1,
                bs=bs*ones(size(s)); % expand bs
            elseif m~=size(s,1) || n~=size(s,2),
                error('bs must be scalar or of same size as s');
            end
            
            x = [];
            y = [];
            for k=1:size(s,2)
                objIndex = bs(1,k);
                [x1,y1] = self.boundary{objIndex}.getXY(s(:,k));
                x(:,k) = x1(:);
                y(:,k) = y1(:);
            end
            if transposed
                x = x.';
                y = y.';
            end
        end
        function self = createBoundaryFunctionCache(self)
            regionToId = self.regionToId;
            self.dlCache = zeros(4,length(self.boundary));
            for k=1:length(self.boundary)
                leftRegionId = regionToId(self.boundary{k}.leftRegion);
                rightRegionId = regionToId(self.boundary{k}.rightRegion);
                self.dlCache(:,k) = [self.boundary{k}.startParam; self.boundary{k}.endParam; leftRegionId; rightRegionId];
            end
        end
    end
    %% Static internal helper methods
    methods(Access = private, Hidden = true, Static)
        function directionDefault = isClockwiseDefault()
            % default direction for closed objects
            directionDefault = true;
        end
        function openDefault = isOpenDefault()
            % default of closed or open for objects defined by a set of
            % lines/arcs
            openDefault = false;
        end
    end
    %% Methods for creating pre-defined geometries
    methods(Static)
        function obj = createCircle(varargin)
            %createCircle Create circle
            %   geometryObject.createCircle(Param1, Value1, Param2, Value2, ...) ... Valid values for
            %   PARAM are:
            %       'name'        : Name of object as char.
            %       'center'      : Coordinates of center point as pointObject.
            %       'radius'      : Radius of circle as scalar double.
            %       'innerRegion' : Name of region inside circle as char.
            %       'outerRegion' : Name of region inside circle as char.
            %   The exteriorRegion property is set to value of 'outerRegion' parameter.
            %
            %   Example: Create a circle named 'circle1' of
            %   radius=0.25 at center=(1,1). The region inside the circle
            %   is called 'regionA'. The region outside the circle is
            %   called 'regionB'.
            %
            %     circle =
            %     geometryObject.createCircle('name','circle1','radius',0.25,'center',pointObject(1,1),'innerRegion','regionA','outerRegion','regionB');
            %
            %   See also pointObject
            import pdetbplus.*;
            ps = parseInputArgsObject.parseInputArgs(varargin{:});
            name = ps.Results.name;
            center = ps.Results.center;
            radius = ps.Results.radius;
            innerRegion = ps.Results.innerRegion;
            outerRegion = ps.Results.outerRegion;
            
            a = cell(0);
            % note arc defined as below is anti-clockwise so we want
            % leftRegion to be the innerRegion
            a{end+1} = arcObject(strcat(name,'Arc1'),'center',center,'radius',radius,'startAngle',0,'endAngle',pi/8);
            a{end}.leftRegion = innerRegion;
            a{end}.rightRegion = outerRegion;
            for k=2:15
                arcName = strcat('Arc',num2str(k));
                a{end+1} = arcObject(strcat(name,arcName),'center',center,'startPoint',a{end}.endPoint(),'rotationAngle',pi/8);
                a{end}.leftRegion = innerRegion;
                a{end}.rightRegion = outerRegion;
            end
            a{end+1} = arcObject(strcat(name,'Arc8'),'center',center,'startPoint',a{end}.endPoint(),'endPoint',a{1}.startPoint());
            a{end}.leftRegion = innerRegion;
            a{end}.rightRegion = outerRegion;
            obj = geometryObject(name,a);
            obj.exteriorRegion = outerRegion;
        end
        function obj = createPolygon(varargin)
            %createPolygon Create polygon
            %   geometryObject.createPolygon(Param1, Value1, Param2, Value2, ...) ... Valid values for
            %   PARAM are:
            %       'name'                 : Name of object as char.
            %       'points'               : Sequential array of coordinates of traversed corner points of polygon as cell array of pointObjects.
            %       'leftRegion'           : Name(s) of region to the left of the direction of traversal of polygon edges as char or cell array of char.
            %       'rightRegion'          : Name(s) of region to the right of the direction of traversal of polygon edges as char or cell array of char.
            %       'leftRegionIsExterior' : Specify if left region is exterior region as bool. Default is true;
            %       'leaveOpen'            : Specify if end and begin points of polygon should be connected as bool. Default is false.
            %
            %   Example: Create a closed polygon named 'poly4' with points:
            %   (0,0), (0,1), (1,1) and (1,-1). The region inside the
            %   polygon is 'regionB' and the region outside the polygon is
            %   'regionA'. The exterior region is 'regionA'.
            %
            %     clockwisePts = ...
            %     {pointObject(0,0),pointObject(0,1),pointObject(1,1),pointObject(1,-1)};
            %     polyg = ...
            %     geometryObject.createPolygon('name','poly4','points',clockwisePts,'leftRegion','regionA','rightRegion','regionB','leftRegionIsExterior',true);
            %
            %   See also pointObject
            % 
            import pdetbplus.*;
            ps = parseInputArgsObject.parseInputArgs(varargin{:});
            name = ps.Results.name;
            pts = ps.Results.points;
            leftRegion = ps.Results.leftRegion;
            rightRegion = ps.Results.rightRegion;
            leftRegionIsExterior = ps.Results.leftRegionIsExterior;
            if isempty(leftRegionIsExterior)
                leftRegionIsExterior = geometryObject.isClockwiseDefault();
            end
            leaveOpen = ps.Results.leaveOpen;
            if isempty(leaveOpen)
                leaveOpen = geometryObject.isOpenDefault();
            end
            % partition with an interior point as specified by
            % interiorPt; Useful for specifying BC edge; creates
            % enumerated dummy interior regions
            interiorPt = ps.Results.interiorPt;
            if ~iscell(leftRegion)
                leftRegionName = leftRegion;
                leftRegion = cell(length(pts),1);
                for k=1:length(leftRegion)
                    leftRegion{k} = leftRegionName;
                end
            end
            if ~iscell(rightRegion)
                rightRegionName = rightRegion;
                rightRegion = cell(length(pts),1);
                for k=1:length(rightRegion)
                    rightRegion{k} = rightRegionName;
                end
            end
            c = cell(0);
            for k=1:length(pts)-1
                nameLine = strcat(name,strcat('Line',num2str(k)));
                c{end+1} = lineObject(nameLine,pts{k},pts{k+1});
                c{end}.leftRegion = leftRegion{k};
                c{end}.rightRegion = rightRegion{k};   
            end
            if ~leaveOpen
                nameLine = strcat(name,strcat('Line',num2str(length(pts))));
                c{end+1} = lineObject(nameLine,pts{end},pts{1});
                c{end}.leftRegion = leftRegion{end};
                c{end}.rightRegion = rightRegion{end};
            end
            if ~isempty(interiorPt)
                c{end+1} = lineObject('somename',pts{1},interiorPt);
                if leftRegionIsExterior
                    c{end}.leftRegion = rightRegion{1};
                    c{end}.rightRegion = rightRegion{end};
                else
                    c{end}.leftRegion = leftRegion{end};
                    c{end}.rightRegion = leftRegion{1};
                end                
                for k=2:length(rightRegion)
                    c{end+1} = lineObject('somename',pts{k},interiorPt);
                    if leftRegionIsExterior
                        c{end}.leftRegion = rightRegion{k};
                        c{end}.rightRegion = rightRegion{k-1};
                    else
                        c{end}.leftRegion = leftRegion{k-1};
                        c{end}.rightRegion = leftRegion{k};
                    end
                end
            end
            obj = geometryObject(name,c);
            if leftRegionIsExterior
                obj.exteriorRegion = leftRegion{1};
            else
                obj.exteriorRegion = rightRegion{1};
            end
        end
        function obj = createSquare(varargin)
            %createSquare Create square
            %   geometryObject.createSquare(Param1, Value1, Param2, Value2, ...) ... Valid values for
            %   PARAM are:
            %       'name'             : Name of object as char.
            %       'edgeStartPoint'   : Coordinates of start point of first edge of square as pointObject
            %       'edgeEndPoint'     : Coordinates of end point of first edge of square as pointObject
            %       'innerRegion'      : Name of region inside of square as char.
            %       'outerRegion'      : Name of region outside of square as char.
            %       'clockwise'        : Specify if traversal of edges of square starting with first edge is clockwise as bool. Default is true;
            %       'leaveOpen'        : Specify if end and begin points of polygon should be connected as bool. Default is false.   
            %   'outerRegion' is always the exterior region.
            %
            %   Example: Create a closed square named 'mysq' with points:
            %   (0,0), (0,1), (1,1) and (1,0). The region inside the
            %   square is 'regionB' and the region outside the square is
            %   'regionA'.
            %
            %     circ = ...
            %     geometryObject.createSquare('name','mysq','edgeStartPoint',pointObject(0,0),'edgeEndPoint',pointObject(0,1),'innerRegion','regionB','outerRegion','regionA','clockwise',true);
            %
            %   See also pointObject
            % 
            import pdetbplus.*;
            ps = parseInputArgsObject.parseInputArgs(varargin{:});
            name = ps.Results.name;
            pt1 = ps.Results.edgeStartPoint;
            pt2 = ps.Results.edgeEndPoint;
            insideRegion = ps.Results.innerRegion;
            outsideRegion = ps.Results.outerRegion;
            clockwise = ps.Results.clockwise;
            leaveOpen = ps.Results.leaveOpen;
            
            if isempty(clockwise)
                clockwise = geometryObject.isClockwiseDefault();
            end
            if isempty(leaveOpen)
                leaveOpen = geometryObject.isOpenDefault();
            end            
            pts = cell(4,0);
            pts{1} = pt1;
            pts{2} = pt2;
            for k = 3:4
                pts{k} = pointObject(pts{k-1}.x-pts{k-2}.x, pts{k-1}.y-pts{k-2}.y) + pts{k-1};
                rotationAngle = pi/2;
                if clockwise
                    rotationAngle = -rotationAngle;
                end
                pts{k} = pts{k}.rotate(pts{k-1},rotationAngle);
            end
            if clockwise
                leftRegion = outsideRegion;
                rightRegion = insideRegion;
                leftRegionIsExterior = true;
            else
                leftRegion = insideRegion;
                rightRegion = outsideRegion;
                leftRegionIsExterior = false;
            end
            obj = geometryObject.createPolygon('name',name,'points',pts,'leftRegion',leftRegion,'rightRegion',rightRegion,...
                'leftRegionIsExterior',leftRegionIsExterior,'leaveOpen',leaveOpen);
        end
        function [objs,isClockwise] = createGeometriesFromImage(varargin)
           %createGeometriesFromImage Create geometries from image
            %   [objs,isClockwise] = geometryObject.createGeometriesFromImage(Param1, Value1, Param2, Value2). Valid values for
            %   PARAM are:
            %       'image'            : Name of image as char.
            %       'minimumPoints'    : Minimum number of points in piecewise linear representation of boundary of geometry as int. Default is 10.
            %   objs                   : Multiple geometries from image as cell array of geometryObject
            %   isClockwise            : Direction of points in piecewise linear representation of boundary of geometry as bool
            %
            %   This method requires a license of **Image Processing Toolbox**.
            %   It is recommended that the method be used along with
            %   geometryObject.setInteriorExteriorRegions(); Prior to calling
            %   geometryObject.setInteriorExteriorRegions(), investigate which
            %   of the image generated geometries are valid.
            %
            %   Example: Create a geometryObject called "img" from one of geometries
            %   generated from 'testImage.png'. Set the number of line boundary
            %   segments on selected geometry to be 15.
            %   The region inside the selected geometry is 'regionB' and the region outside the selected geometry is
            %   'regionA'.
            %
            %     [objs,isClockwise] = ...
            %     geometryObject.createGeometriesFromImage('image','testImage.png','minimumPoints',15);
            %     % objs{3} corresponds to an expected geometry. use objs{3}.plot() to visually inspect geometry.
            %     img = objs{3}.setInteriorExteriorRegions('regionB','regionA',isClockwise);
            %
            %   See also geometryObject/setInteriorExteriorRegions
            %             
            import pdetbplus.*;
            ps = parseInputArgsObject.parseInputArgs(varargin{:});
            imageName = ps.Results.image;
            minPts = ps.Results.minimumPoints;
            if isempty(minPts)
                minPts = 10;
            end
            im = imread(imageName);
            bw = im2bw(im);
            boundaries = bwboundaries(bw);
            numBoundaries = length(boundaries);
            objs = cell(numBoundaries,1);
            for k=1:numBoundaries
                b = boundaries{k};
                numPts = size(b,1);
                if numPts < minPts
                    continue;
                end
                b = double(b);
                quarterPoints = floor(numPts/4);
                pwl{1} = piecewiseLineObject(strcat(imageName,num2str(4*k)),b(1:quarterPoints,2),-b(1:quarterPoints,1));
                pwl{1}.leftRegion = 'dummyRegionA';
                pwl{1}.rightRegion = 'dummyRegionB';
                pwl{2} = piecewiseLineObject(strcat(imageName,num2str(4*k+1)),b(quarterPoints:2*quarterPoints,2),-b(quarterPoints:2*quarterPoints,1));
                pwl{2}.leftRegion = 'dummyRegionA';
                pwl{2}.rightRegion = 'dummyRegionB';
                pwl{3} = piecewiseLineObject(strcat(imageName,num2str(4*k+2)),b(2*quarterPoints:3*quarterPoints,2),-b(2*quarterPoints:3*quarterPoints,1));
                pwl{3}.leftRegion = 'dummyRegionA';
                pwl{3}.rightRegion = 'dummyRegionB';
                pwl{4} = piecewiseLineObject(strcat(imageName,num2str(4*k+3)),[b(3*quarterPoints:end,2);b(1,2)],-[b(3*quarterPoints:end,1);b(1,1)]);
                pwl{4}.leftRegion = 'dummyRegionA';
                pwl{4}.rightRegion = 'dummyRegionB';
                objs{k} = geometryObject(strcat(imageName,num2str(k)),pwl);
                objs{k}.exteriorRegion = 'dummyRegionA';
            end
            isClockwise = true; % AKA clockwise
        end
    end
end
