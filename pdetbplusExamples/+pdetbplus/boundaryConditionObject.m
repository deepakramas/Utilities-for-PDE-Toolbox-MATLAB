classdef boundaryConditionObject < pdetbplus.formulation
    %boundaryConditionObject Create boundary condition (BC) object
    %   Use boundaryConditionObject to specify BCs, one at a time, and pass BC function handle
    %   to PDE Toolbox functions assempde, adaptmesh, parabolic, hyperbolic.
    %
    %   boundaryConditionObject methods:
    %
    %       add         - Add Boundary condition, 1 boundary at a time
    %       getMatrices - Get Q,G,H,R matrices as described in PDE Toolbox function assemb
    %       getGlobalHR - Get global H and R matrix contributions
    %
    %   boundaryConditionObject properties:
    %
    %       bcFunction  - BC function with identical purpose and signature as PDE Toolbox function pdebound
    %
    %   Example: create boundaryConditionObject for geometry, g1 and problem output
    %   dimensions = 2.
    %
    %       bc = boundaryConditionObject(g1,2);
    %
    %   Example: add Dirichlet Boundary condition of [0,0,1] to boundaryConditionObject, bc for
    %   boundary named 'A' in 3 output dimensional problem.
    %
    %       bc.add('name','A','dirichlet',[0,0,1]);
    %
    %   Example: specify Neumann BC (= t in 3rd dimension, 0 otherwise) in 3 output dimensional problem
    %   on boundary between regions 'A' and 'B' and add bc to boundaryConditionObject, bc.
    %
    %       function [hval,rval,qval,gval] = foo(x,y,u,t)
    %           N = 3; % #output dimensions
    %           rval = []; hval = [];
    %           qval = sparse(N,N); gval = [0;0;t];
    %       end
    %       bc.add('innerRegion','A','outerRegion','B','xyutFunction',@foo);
    %
    %   Example: pass boundary condition function from boundaryConditionObject bc to assempde.
    %
    %       assempde(bc.bcFunction,...);
    %
    %   See also geometryObject pdebound assemb assempde adaptmesh parabolic hyperbolic coeffsObject
    %
    % 
    properties (Access = protected, Hidden = true)
        bc = cell(0);
        parser = pdetbplus.parseInputArgsObject; 
    end
    properties (Dependent = true,SetAccess = private)
        %bcFunction BC function with identical purpose and signature as PDE Toolbox function pdebound
        %
        %   Example: pass boundary condition function from boundaryConditionObject bc to assempde.
        %       assempde(bc.bcFunction,...);
        %
        %   See also pdebound.
        bcFunction;
    end
    methods
        function self = boundaryConditionObject(geometry, dimension)
            %boundaryConditionObject Create BC object
            %   boundaryConditionObject(geometry, dimension) returns a boundaryConditionObject
            %   instance for input: geometryObject, geometry and dimension = #output Dimension of
            %   problem.
            %
            %   Example: create boundaryConditionObject for geometry, g1 and problem output
            %   dimensions = 2.
            %       bc = boundaryConditionObject(g1,2);
            %
            %   See also geometryObject coeffsObect.
            self = self@pdetbplus.formulation(geometry, dimension); 
        end       
        function add(self,varargin)
            %add Add Boundary condition, 1 boundary at a time
            %   add(PARAM1,VALUE1,PARAM2,VALUE2,...) adds a boundary condition to
            %   boundaryConditionObject instance. PARAM can be:
            %
            %       'name' : Name of boundary to which BC applies as char.
            %
            %       'innerRegion','outerRegion' : Names of inner and outer regions of boundary as
            %       char to which BC applies. The boundary edge normal points from 'innerRegion'
            %       outwards to 'outerRegion'. 'name' takes precedence over 'innerRegion' and
            %       'outerRegion'.
            %
            %       'dirichlet', 'neumann'      : Value of Dirichlet and Neumann BC as an array of
            %       size #Output dimensions X 1.
            %          
            %       'xyutFunction'              : Value of boundary function as function handle that
            %       is of the form foo(x,y,u,time) where x,y,u are boundary edge point values and u
            %       is of length #Output dimensions. x,y are position coordinates and u is the
            %       output that is being solved for. The format for xyutFunction is:
            %           function [hval,rval,qval,gval] = fooBcondFunction(x,y,u,t)
            %               rval = N x 1 vector or [];
            %               hval = N x N matrix or [];
            %               qval = N x N matrix or [];
            %               gval = N x 1 vector or [];
            %           end
            %       See interpretation of hval,rval,qval,gval in PDE Toolbox function pdebound. Note
            %       that hval,rval,qval and gval are atomic valued i.e. for a single point on the
            %       boundary. hval,rval,qval,gval allow for coupling between variables corresponding
            %       to all output dimensions at that node only i.e. if u1(x1,y1),u2(x1,y1) are two
            %       output variables and N = 2 => at a node(x1,y1), u1(x1,y1) and u2(x2,y2) can be coupled
            %       using Dirichlet or Generalized Neumann expressions.
            %
            %       Example: add Dirichlet Boundary condition of [0,0,1] to boundaryConditionObject,
            %       bc for boundary named 'A' in 3 output dimensional problem.
            %
            %           bc.add('name','A','dirichlet',[0,0,1]);
            %
            %       Example: add Neumann Boundary condition of [0,1,1] to boundaryConditionObject,
            %       bc for boundary between regions 'A' and 'B' in 3 output dimensional problem.The
            %       boundary normal points from 'A' to 'B'.
            %
            %           bc.add('innerRegion','A','outerRegion',B,'neumann',[0,0,1]);
            %
            %       Example: specify Neumann BC (= t in 3rd dimension, 0 otherwise) in 3 output dimensional problem
            %       on boundary between regions 'A' and 'B' and add bc to boundaryConditionObject, bc.
            %
            %           function [hval,rval,qval,gval] = foo(x,y,u,t)
            %               N = 3; % #output dimensions
            %               rval = []; hval = [];
            %               qval = sparse(N,N); gval = [0;0;t];
            %           end
            %           bc.add('innerRegion','A','outerRegion','B','xyutFunction',@foo);
            %
            %   See also pdebound bcFunction getGlobalHR getMatrices assemb.
            ps = self.parser.parseInputArgs(varargin{:});
            self.bc{end+1}.xyutFunction = ps.Results.xyutFunction;
            self.bc{end}.outerRegion = ps.Results.outerRegion;
            self.bc{end}.innerRegion = ps.Results.innerRegion;
            self.bc{end}.name = ps.Results.name;

            if isempty(self.bc{end}.xyutFunction)
                if ~isempty(ps.Results.dirichlet)
                    self.bc{end}.xyutFunction = @(x,y,u,time) self.dirichlet(ps.Results.dirichlet,self.dimension);
                elseif ~isempty(ps.Results.neumann)
                    self.bc{end}.xyutFunction = @(x,y,u,time) self.neumann(ps.Results.neumann,self.dimension);
                end
            end
            for bcix=1:length(self.bc)
                self.bc{bcix}.boundaryMatchIndex = [];
                if ~isempty(self.bc{bcix}.name) && ~isempty(self.bc{bcix}.xyutFunction)
                    for bcix2 = 1:length(self.geometry.boundary)
                        if strcmp(self.geometry.boundary{bcix2}.name,self.bc{bcix}.name)
                            self.bc{bcix}.boundaryMatchIndex = bcix2;
                        end
                    end
                end
            end
        end
        function bcf = get.bcFunction(self)            
            bcf = @self.bcFunction_impl;
        end
        function [Q,G,H,R] = getMatrices(self,varargin)
            %getMatrices Get Q,G,H,R matrices as described in PDE Toolbox function assemb
            %   [Q,G,H,R] = getMatrices(PARAM1,VALUE1,PARAM2,VALUE2,...) where PARAM can be:
            %       'u'   : Output variable of PDE as an array (#Dimensions X #Nodes) long.
            %       'time': Value of time as double.
            %       'p','e','t': Values of externally specified mesh of structure as described in initmesh.
            %
            %   Example: Get Q,G,H,R matrices corresponding to u = u0 for boundaryConditionObject, bc.
            %       [Q,G,H,R] = bc.getMatrices('u',u0);
            %
            %   See also assemb initmesh meshObject.
            ps = self.parser.parseInputArgs(varargin{:});
            u = ps.Results.u;
            time = ps.Results.time;
            p = ps.Results.p;
            e = ps.Results.e;
            if isempty(p)
                p = self.geometry.mesh.p;
            end
            if isempty(e)
                e = self.geometry.mesh.e;
            end

            if isempty(u) && isempty(time)
                [Q,G,H,R] = assemb(self.bcFunction,p,e);
            elseif isempty(time)
                [Q,G,H,R] = assemb(self.bcFunction,p,e,u);
            elseif isempty(u)
                [Q,G,H,R] = assemb(self.bcFunction,p,e,time);
            else
                [Q,G,H,R] = assemb(self.bcFunction,p,e,u,time);
            end
        end       
        function [Hadd,Radd] = getGlobalHR(self,varargin)
            %getGlobalHR Get global H and R matrix contributions
            %   [Hadd,Radd] = getGlobalHR(PARAM1,VALUE1,PARAM2,VALUE2,...) returns H,R matrix (see
            %   assemb in PDE Toolbox) contributions corresponding to 'global' BC such as periodic,
            %   tie constraints etc, where PARAM can be:
            %
            %       'name'   : Name of a boundary or a region as char. If a boundary and a region have
            %       the same name, the boundary is selected. Specifiying 'name' identifies all the nodes
            %       on the boundary/region and further operations are performed as described below.
            %
            %       'nodes'  : External list of mesh node indices as array of int. The nodes
            %       corresponding to 'name' take precedence if 'name' is specified.
            %
            %       'dimension' : Dimension of output variable where the global BC is applicable as int.
            %
            %       'value'  : Value of the boundary condition as double.
            %
            %       'type'   : Has value of 'tie' or 'point' or 'periodic'. 'tie' forces the prescribed nodes to
            %       have the same value. 'point' constrains the prescribed nodes to have the same
            %       'value'. 'periodic' imples prescribed nodes have the same output value as
            %       another boundary described by 'periodicX' and 'periodicY'.
            %
            %       'periodicX','periodicY' : Curve coordinates as char. Applicable only if 'type'
            %       is 'periodic'. Curve coordinates are expressed in terms of coordinates of points
            %       on the reference nodes - 'x' and 'y'. The curve should be inside or on the
            %       geometry.
            %
            %       'showNodes' : Plot nodes where BCs are applied. Default is false.
            %
            %   H,R matrix contributions from various constraints can be combined as H = [H1;H2;...]
            %   and R = [R1;R2;...] and used in assempde.
            %       Example: Tie the 3rd dimension of output for boundary 'A1' and return H,R matrix
            %       contributions for boundaryConditionObject, bc.
            %
            %           [H1,R1] = bc.getGlobalHR('name','A1','dimension',3,'type','tie');
            %
            %       Example: Constrain the 1st dimension of output for region 'B' to have the same
            %       value, 2.5 for boundaryConditionObject, bc.
            %
            %           [H1,R1] = bc.getGlobalHR('name','B','dimension',1,'type','point','value',2.5);
            %
            %       Example: Constrain the list of nodes n1 to have the same value along output
            %       dimension 1 for boundaryConditionsObject, bc.
            %
            %           [H1,R1] = bc.getGlobalHR('nodes',n1,'dimension',3,'type','tie');
            %
            %       Example: Apply periodic bounditions on boundary 'B' and boundary 'x+3' where 'x'
            %       is a point on boundary 'B'. The periodic BCs are applied along output dimension 1 for
            %       boundaryConditionsObject, bc.
            %
            %           [H1,R1] =
            %           bc.getGlobalHR('name','B','dimension',1,'type','periodic','periodicX','x+3',periodicY,'0');
            %
            %   See also assemb assempde.
            ps = self.parser.parseInputArgs(varargin{:});
            name = ps.Results.name;
            type = ps.Results.type;
            value = ps.Results.value;
            if isempty(value)
                value = 0;
            end
            bNodes = ps.Results.nodes;
            N = ps.Results.dimension;
            if isempty(N)
                error 'dimension must be specified';
            end

            numNodes = size(self.geometry.mesh.p,2);
            if ~isempty(name)
                bNodes = self.geometry.getBoundaryNodes('name',name);
                if isempty(bNodes)
                    bNodes = self.geometry.getRegionNodes('region',name);
                end
            end
            switch type
                case 'tie'                    
                    Hadd = sparse(length(bNodes)-1,self.dimension*numNodes);
                    Radd = ones(size(Hadd,1),1)*value;
                    for kk = 2:length(bNodes)
                        Hadd(kk-1,(N-1)*numNodes + bNodes(kk-1)) = 1;
                        Hadd(kk-1,(N-1)*numNodes + bNodes(kk)) = -1;
                    end
                    if ps.Results.showNodes
                        x = self.geometry.mesh.p(1,bNodes)';
                        y = self.geometry.mesh.p(2,bNodes)';
                        plot(x,y,'*r');
                    end
                case 'point'
                    Hadd = sparse(1,self.dimension*numNodes);
                    Radd = zeros(1,1);
                    Hadd(1,(N-1)*numNodes + bNodes) = 1;
                    Radd(1) = value;
                    if ps.Results.showNodes
                        x = self.geometry.mesh.p(1,bNodes)';
                        y = self.geometry.mesh.p(2,bNodes)';
                        plot(x,y,'*r');
                    end
                case 'gradient'
                    import pdetbplus.xyFunctionClass;                    
                    [dudxContrib,dudyContrib] = xyFunctionClass([],self.geometry.mesh).gradientContributions();
                    Hadd = sparse(size(dudxContrib,1)+size(dudyContrib,1), self.dimension*numNodes);
                    Hadd(:,((N-1)*numNodes + 1):N*numNodes) = [dudxContrib;dudyContrib];
                    Radd = zeros(size(Hadd,1),1);
                    Radd = value;
                case 'periodic'
                    periodicX = ps.Results.periodicX;
                    periodicY = ps.Results.periodicY;
                    x = self.geometry.mesh.p(1,bNodes)';
                    y = self.geometry.mesh.p(2,bNodes)';
                    if isempty(periodicX) || isempty(periodicY)
                        error 'for periodic type, periodicXShift AND periodicYShift expressions must be specified';
                    end
                    newX = eval(periodicX);
                    newY = eval(periodicY);
                    if ps.Results.showNodes
                        plot(x,y,'*b');plot(newX,newY,'*r');
                        quiver(x,y,newX-x,newY-y,0,'--k','filled');
                    end
                    import pdetbplus.xyFunctionClass;
                    nodesAndContributions = xyFunctionClass([],self.geometry.mesh).elementIndexFunction(newX,newY);
                    if find(isnan(nodesAndContributions(:,1)))
                        error 'At least some noes of curve mapped by periodicX and periodicY are not in geometry';
                    end
                    Hadd = sparse(length(bNodes),self.dimension*numNodes);
                    Radd = ones(size(Hadd,1),1)*value;
                    for kk = 1:length(bNodes)
                        Hadd(kk,(N-1)*numNodes + bNodes(kk)) = 1;
                        node1 = nodesAndContributions(kk,1);
                        node2 = nodesAndContributions(kk,2);
                        node3 = nodesAndContributions(kk,3);
                        Hadd(kk,(N-1)*numNodes + node1) = -nodesAndContributions(kk,4);
                        Hadd(kk,(N-1)*numNodes + node2) = -nodesAndContributions(kk,5);
                        Hadd(kk,(N-1)*numNodes + node3) = -nodesAndContributions(kk,6);
                    end
            end
        end
    end
    methods (Static,Access = private,Hidden = true)
        function [hval,rval,qval,gval] = dirichlet(rval,N)
            % Dirichlet condition on the boundary
            hval = speye(N);
            qval = [];
            gval = [];
        end
        function [hval,rval,qval,gval] = neumann(gval,N)
            % Neumann condition on the boundary
            rval = [];
            hval = [];
            qval = sparse(N,N);
        end
    end
    methods (Access = private, Hidden = true)
        function [q,g,h,r] = bcFunction_impl(self,p,e,u,time)
            ne = size(e,2); % number of edges
            N = self.dimension;
            q = zeros(N^2,ne);
            g = zeros(N,ne);
            h = zeros(N^2,2*ne);
            r = zeros(N,2*ne);
            for ix=1:ne
                xyutFunction = self.getXyutFunctionRegions(e(self.geometry.mesh.edgeLeftRegionIndex,ix),e(self.geometry.mesh.edgeRightRegionIndex,ix));
                if isempty(xyutFunction)
                    p1Index = e(self.geometry.mesh.edgeFirstPointIndex,ix); % note e is just the boundary edges
                    p2Index = e(self.geometry.mesh.edgeSecondPointIndex,ix); % note e is just the boundary edges
                    originalp1 = pdetbplus.pointObject(p(1,p1Index),p(2,p1Index));
                    originalp2 = pdetbplus.pointObject(p(1,p2Index),p(2,p2Index));
                    xyutFunction = self.getXyutFunctionName(originalp1,originalp2);
                end
                if ~isempty(xyutFunction)
                    p1 = e(self.geometry.mesh.edgeFirstPointIndex,ix);
                    p2 = e(self.geometry.mesh.edgeSecondPointIndex,ix);
                    if isempty(u)
                        up1 = [];
                        up2 = [];
                    else
                        up1 = u(p1,:);
                        up2 = u(p2,:);
                    end
                    if nargin(xyutFunction) == 4
                        [h1,r1,~,~] = xyutFunction(p(1,p1),p(2,p1),up1,time);
                        [h2,r2,~,~] = xyutFunction(p(1,p2),p(2,p2),up2,time);
                        [~,~,q1,g1] = xyutFunction(1/2*(p(1,p2)+p(1,p1)),1/2*(p(2,p2)+p(2,p1)),up2,time);
                    elseif nargin(xyutFunction) == 5
                        [h1,r1,~,~] = xyutFunction(p(1,p1),p(2,p1),up1,time,ix);
                        [h2,r2,~,~] = xyutFunction(p(1,p2),p(2,p2),up2,time,ix);
                        [~,~,q1,g1] = xyutFunction(1/2*(p(1,p2)+p(1,p1)),1/2*(p(2,p2)+p(2,p1)),up2,time,ix);                        
                    end
                    if ~isempty(h1) && ~isempty(h2) && ~isempty(r1) && ~isempty(r2)
                        h(:,ix) = h1(:);
                        h(:,ix+ne) = h2(:);
                        r(:,ix) = r1(:);
                        r(:,ix+ne) = r2(:);
                    end
                    if ~isempty(q1) && ~isempty(g1)
                        q(:,ix) = q1(:);
                        g(:,ix) = g1(:);
                    end
                end
            end
        end
        function xyutFunction = getXyutFunctionRegions(self,region1Id,region2Id)
            regionToIds = self.geometry.regionToId;
            xyutFunction = [];
            for bcix=1:length(self.bc)
                if ~isempty(self.bc{bcix}.outerRegion) && ~isempty(self.bc{bcix}.innerRegion)
                    bcOuterRegionId = regionToIds(self.bc{bcix}.outerRegion);
                    bcInnerRegionId = regionToIds(self.bc{bcix}.innerRegion);
                    if ((region1Id == bcOuterRegionId) && (region2Id == bcInnerRegionId)) || ((region1Id == bcInnerRegionId) && (region2Id == bcOuterRegionId))
                        xyutFunction = self.bc{bcix}.xyutFunction;
                        break;
                    end
                end
            end
        end
        function xyutFunction = getXyutFunctionName(self,p1,p2)
            xyutFunction = [];
            for bcix=1:length(self.bc)
                if ~isempty(self.bc{bcix}.boundaryMatchIndex)
                    boundaryMatchIndex = self.bc{bcix}.boundaryMatchIndex;
                    boundary = self.geometry.boundary{boundaryMatchIndex};
                    if ismethod(boundary,'isOn') && boundary.isOn(p1) && boundary.isOn(p2)
                        xyutFunction = self.bc{bcix}.xyutFunction;
                        break;
                    end
                end
            end
        end
    end
end