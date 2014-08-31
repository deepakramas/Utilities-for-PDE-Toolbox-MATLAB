classdef coeffsObject < pdetbplus.formulation
    %coeffsObject Represent coefficients (c,f,a,d) for geometry
    %   Use coeffsObject to specify coefficients and pass coefficients function handle
    %   to PDE Toolbox functions assempde, adaptmesh, parabolic, hyperbolic.
    %
    %   coeffsObject methods:
    %
    %       add         - Add coefficient per region of geometry
    %       getMatrices - Get area integral matrices as described in PDE Toolbox function assema
    %       getJacobian - Get Jacobian of area integral matrices
    %
    %   coeffsObject properties:
    %
    %       fFunction   - Function handle for "f" coefficient as described in PDE Toolbox doc
    %       cFunction   - Function handle for "c" coefficient as described in PDE Toolbox doc
    %       aFunction   - Function handle for "a" coefficient as described in PDE Toolbox doc
    %       dFunction   - Function handle for "d" coefficient as described in PDE Toolbox doc
    %
    %   Example: create coeffsObject for geometry, g1 and problem output
    %   dimensions = 2.
    %       coeffs = coeffsObject(g1,2);
    %
    %   Example: for a 3 output dimension problem, specify "f" as [0,0,1] for region 'A' in
    %   coeffsObject, coeffs.
    %           function fvalue = fCoeff(x,y,u,ux,uy,time)
    %               fvalue = [0,0,1];
    %           end
    %           coeffs.add('region','A','fiFunction',@fCoeff);
    %
    %   Example: for a scalar output problem, specify "c" for PDE term: [d/dx, d/dy]*[1,
    %   -3.1;0, 4.7]*[du/dx; du/dy] for region 'A' in coeffsObject, coeffs.
    %           function Cvalue = cCoeff(x,y,u,ux,uy,time)
    %               Cvalue = [1,-3.1;0,4.7];
    %           end
    %           coeffs.add('region','A','cijFunction',@cCoeff);
    %
    %   Example: for a 2 output dimension problem, specify "a" for PDE term: [1 0; -1
    %   3.3]*[u1;u2] for region 'A' in coeffsObject, coeffs.
    %           function avalue = aCoeff(x,y,u,ux,uy,time)
    %               avalue = [1,0; -1,3.3];
    %           end
    %           coeffs.add('region','A','aijFunction',@aCoeff);
    %
    %   Example: for a scalar output problem, specify "f" as 2 for region 'A' in
    %   coeffsObject, coeffs.
    %           coeffs.add('region','A','fConstantValue',2);
    %
    %   Example: for a 2 output dimension problem, specify "c" for PDE term [d/dx,d/dy]*[1 0
    %   0 0;0 1 0 0;0 0 1 0;0 0 0 1]*[du1dx;du1dy;du2dx;du2dy] for region 'A' in
    %   coeffsObject, coeffs.
    %           coeffs.add('region','A','cConstantValue',1);
    %
    %   Example: symbolically prescribe scalar Poisson's equation for region, 'A' in
    %   coeffsObject, coeffs.
    %           [equations,variables] = function symbolicPoisson()
    %               syms u ddx ddy dudx dudy;
    %               rho = 1.0;
    %               equations = -[ddx ddy]*[rho 0;0 rho]*[dudx;dudy] + [1 2]; % the PDE
    %               variables = u;
    %           end
    %           coeffs.add('region','A','symbolicEquationFunction',@symbolicPoisson);
    %
    %   See also geometryObject http://www.mathworks.com/help/pde/ug/multidimensional-coefficients.html http://www.mathworks.com/help/pde/pde-coefficients.html assema assempde adaptmesh parabolic hyperbolic boundaryConditionObject
    %
    % 
    properties (Access = private)
        coeff = cell(0);
        triangleIndividualCoordinatesForm = false;
        displayCoefficients = false;
        parser = pdetbplus.parseInputArgsObject;
    end
    properties (Dependent = true,SetAccess = private)
        %fFunction Function handle for "f" coefficient as described in PDE Toolbox doc
        %   fFunction can be passed to PDE Toolbox functions assempde, adaptmesh, parabolic, hyperbolic
        %
        %   Example: pass fFunction from coeffObject, coeffs, to PDE Toolbox function assempde
        %       assempde(...,coeffs.fFunction,...);
        %
        %   See also assempde adaptmesh parabolic hyperbolic
        fFunction;
        %cFunction Function handle for "c" coefficient as described in PDE Toolbox doc
        %   cFunction can be passed to PDE Toolbox functions assempde, adaptmesh, parabolic, hyperbolic
        %
        %   Example: pass cFunction from coeffObject, coeffs, to PDE Toolbox function assempde
        %       assempde(...,coeffs.cFunction,...);
        %
        %   See also assempde adaptmesh parabolic hyperbolic        
        cFunction;
        %aFunction Function handle for "a" coefficient as described in PDE Toolbox doc
        %   aFunction can be passed to PDE Toolbox functions assempde, adaptmesh, parabolic, hyperbolic
        %
        %   Example: pass aFunction from coeffObject, coeffs, to PDE Toolbox function assempde
        %       assempde(...,coeffs.aFunction,...);
        %
        %   See also assempde adaptmesh parabolic hyperbolic         
        aFunction;
        %dFunction Function handle for "d" coefficient as described in PDE Toolbox doc
        %   dFunction can be passed to PDE Toolbox functions assempde, adaptmesh, parabolic, hyperbolic
        %
        %   Example: pass dFunction from coeffObject, coeffs, to PDE Toolbox function assempde
        %       assempde(...,coeffs.dFunction,...);
        %
        %   See also assempde adaptmesh parabolic hyperbolic         
        dFunction;
    end
    methods
        function self = coeffsObject(geometry, dimension)
            %coeffsObject Create PDE coefficients object
            %   coeffsObject(geometry, dimension) returns a coeffsObject
            %   instance for input: geometryObject, geometry and dimension = #output Dimension of
            %   problem.
            %
            %   Example: create coeffsObject for geometry, g1 and problem output
            %   dimensions = 2.
            %       coeffs = coeffsObject(g1,2);
            %
            %   See also geometryObject boundaryConditionObject.
            self = self@pdetbplus.formulation(geometry, dimension); 
        end
        function add(self,varargin)
            %add Add coefficient per region of geometry
            %   add(PARAM1,VALUE1,PARAM2,VALUE2,...) adds coefficients for a region to
            %   coeffsObject instance. PARAM can be:
            %
            %       'region' : Name of region to which coefficient applies as char.
            %
            %       'fiFunction' : Function handle for setting the "f" coefficient value at element centroid. The format
            %       for fiFunction is -
            %           function fi = fCoeff(x,y,u,ux,uy,time) % receive x,y,y,ux,uy at centroids
            %               N = #output dimensions
            %               % x,y are the position coordinates. u is the solution, ux and uy are the
            %               % direction derivatives of u.
            %               % fi is an N length vector; nth row corresponds to entry for nth
            %               % output equation
            %               fi = N X 1 vector; % Can be sparse.
            %           end
            %
            %       'cijFunction' : Function handle for setting the "c" coefficient value at element centroid. The format
            %       for cijFunction is -
            %           function C = cCoeff(x,y,u,ux,uy,time) % receive x,y,y,ux,uy at centroids
            %               % C is an 2N x 2N matrix. C is an N-by-N matrix of 2-by-2 blocks as per the PDE
            %               % Toolbox documentation ("logical" 2N X 2N C matrix -
            %               % http://www.mathworks.com/help/pde/ug/c.html)
            %               C = 2N X 2N matrix; % Can be sparse
            %           end
            %
            %       'aijFunction' OR 'dijFunction' : Function handle for setting the "a" or "d"
            %       coefficient at element centroid. The format for the aijFunction/dijFunction is -
            %           function aij = aCoeff(x,y,u,ux,uy,time) % receive x,y,y,ux,uy at centroids
            %               % aij is an N x N matrix; (m,n) entry in aij corresponds to the
            %               % contribution of the nth output dimension variable to mth output
            %               % equation
            %               aij = N X N matrix. % Can be sparse
            %           end
            %
            %       'fConstantValue' : Scalar or vector value for the "f" coefficient. In the scalar
            %       case, the same value is applied for all output dimensions. The vector value is
            %       of size #output dimensions x 1.
            %
            %       'cConstantValue', 'aConstantValue', 'dConstantValue' : Scalar or vector value
            %       for the diagonal of the "c"/"a"/"d" coefficient matrix. In the scalar case, the
            %       value is set to be same for all diagonal terms e.g. d/dx (.) d/dx and d/dy (.)
            %       d/dy for 'cConstantValue' for all output dimensions where (.) denotes the
            %       coefficient. The vector value is the diagonal of the coefficient matrix and is
            %       #output dimensions X 1 long for 'aConstantValue' and 'dConstantValue' and
            %       2*#output dimensions X 1 long for 'cConstantValue';
            %
            %       'symbolicEquationFunction' : Function handle that describes equations
            %       symbolically using the Symbolic Math Toolbox. A license of Symbolic Math Toolbox
            %       is required. Format for is -
            %           function [equations,variables] = fooSymbolicEquations
            %               equations = array of symbolic equations representing PDE;
            %               variables = variables of PDE;
            %               % Reserved keywords: x,y,time,ddx,ddy,d[variable]dx,d[variable]dy
            %               % Specify the Div operator as [ddx ddy]
            %           end
            %
            %   Example: for a 3 output dimension problem, specify "f" as [0,0,1] for region 'A' in
            %   coeffsObject, coeffs.
            %           function fvalue = fCoeff(x,y,u,ux,uy,time)
            %               fvalue = [0,0,1];
            %           end
            %           coeffs.add('region','A','fiFunction',@fCoeff);
            %
            %   Example: for a scalar output problem, specify "c" for PDE term: [d/dx, d/dy]*[1,
            %   -3.1;0, 4.7]*[du/dx; du/dy] for region 'A' in coeffsObject, coeffs.
            %           function Cvalue = cCoeff(x,y,u,ux,uy,time)
            %               Cvalue = [1,-3.1;0,4.7];
            %           end
            %           coeffs.add('region','A','cijFunction',@cCoeff);
            %
            %   Example: for a 2 output dimension problem, specify "a" for PDE term: [1 0; -1
            %   3.3]*[u1;u2] for region 'A' in coeffsObject, coeffs.
            %           function avalue = aCoeff(x,y,u,ux,uy,time)
            %               avalue = [1,0; -1,3.3];
            %           end
            %           coeffs.add('region','A','aijFunction',@aCoeff);
            %
            %   Example: for a scalar output problem, specify "f" as 2 for region 'A' in
            %   coeffsObject, coeffs.
            %           coeffs.add('region','A','fConstantValue',2);
            %
            %   Example: for a 2 output dimension problem, specify "c" for PDE term [d/dx,d/dy]*[1 0
            %   0 0;0 1 0 0;0 0 1 0;0 0 0 1]*[du1dx;du1dy;du2dx;du2dy] for region 'A' in
            %   coeffsObject, coeffs.
            %           coeffs.add('region','A','cConstantValue',1);
            %
            %   Example: symbolically prescribe scalar Poisson's equation for region, 'A' in
            %   coeffsObject, coeffs.
            %           [equations,variables] = function symbolicPoisson()
            %               syms u ddx ddy dudx dudy;
            %               rho = 1.0;
            %               equations = -[ddx ddy]*[rho 0;0 rho]*[dudx;dudy] + [1 2]; % the PDE
            %               variables = u;
            %           end
            %           coeffs.add('region','A','symbolicEquationFunction',@symbolicPoisson);
            %
            %   See also http://www.mathworks.com/help/pde/ug/c.html http://www.mathworks.com/help/pde/ug/multidimensional-coefficients.html
            ps = self.parser.parseInputArgs(varargin{:});
            self.coeff{end+1}.region = ps.Results.region;
            if ~isempty(ps.Results.fConstantValue)                
                fConstantValue = ps.Results.fConstantValue;
                fiConstantFunction = @(xc,yc,uc,uxc,uyc,time) self.fiConstantFunctionHelper(fConstantValue,self.dimension);
                self.coeff{end}.fiFunction = fiConstantFunction;
            else
                self.coeff{end}.fiFunction = ps.Results.fiFunction;
            end
            if ~isempty(ps.Results.cConstantValue)
                cConstantValue = ps.Results.cConstantValue;
                cijConstantDiagFunction = @(xc,yc,uc,uxc,uyc,time) self.cijConstantDiagFunctionHelper(cConstantValue,self.dimension);
                self.coeff{end}.cijFunction = cijConstantDiagFunction;
            else
               self.coeff{end}.cijFunction = ps.Results.cijFunction; 
            end
            if ~isempty(ps.Results.aConstantValue)
                aConstantValue = ps.Results.aConstantValue;
                aijConstantDiagFunction = @(xc,yc,uc,uxc,uyc,time) self.aijConstantDiagFunctionHelper(aConstantValue,self.dimension);
                self.coeff{end}.aijFunction = aijConstantDiagFunction;
            else
               self.coeff{end}.aijFunction = ps.Results.aijFunction;
            end
            if ~isempty(ps.Results.dConstantValue)
                dConstantValue = ps.Results.dConstantValue;
                dijConstantDiagFunction = @(xc,yc,uc,uxc,uyc,time) self.aijConstantDiagFunctionHelper(dConstantValue,self.dimension);
                self.coeff{end}.dijFunction = dijConstantDiagFunction;
            else
               self.coeff{end}.dijFunction = ps.Results.dijFunction;
            end            
            
            self.coeff{end}.symbolicEquationFunction = ps.Results.symbolicEquationFunction; 
            self.coeff{end}.jcFunction = [];
            self.coeff{end}.jaFunction = [];
            self.coeff{end}.jdFunction = [];
            self.coeff{end}.f = cell(self.dimension,1);
            
            if ~isempty(self.coeff{end}.symbolicEquationFunction)
                [equations,variables] = self.coeff{end}.symbolicEquationFunction();
                [mfC,mfA,mfF,mfD,jfC,jfA,jfD] = self.generateMatrixFunctionsFromSymbolic('index',length(self.coeff),'equations',equations,'variables',variables,'displayCoefficients',self.displayCoefficients);
                if isempty(self.coeff{end}.cijFunction) && ~isempty(mfC)
                    tempF = @(x,y,u,ux,uy,time) self.coeffFromSymbolicEquation(mfC,2,self.dimension,x,y,u,ux,uy,time);
                    self.coeff{end}.cijFunction = tempF;
                    tempF = @(x,y,u,ux,uy,time) self.coeffFromSymbolicEquation(jfC,2,self.dimension,x,y,u,ux,uy,time);
                    self.coeff{end}.jcFunction = tempF;
                end
                if isempty(self.coeff{end}.aijFunction) && ~isempty(mfA)
                    tempF = @(x,y,u,ux,uy,time) self.coeffFromSymbolicEquation(mfA,1,self.dimension,x,y,u,ux,uy,time);
                    self.coeff{end}.aijFunction = tempF;
                    tempF = @(x,y,u,ux,uy,time) self.coeffFromSymbolicEquation(jfA,1,self.dimension,x,y,u,ux,uy,time);
                    self.coeff{end}.jaFunction = tempF;                    
                end                
                if isempty(self.coeff{end}.dijFunction) && ~isempty(mfD)
                    tempF = @(x,y,u,ux,uy,time) self.coeffFromSymbolicEquation(mfD,1,self.dimension,x,y,u,ux,uy,time);
                    self.coeff{end}.dijFunction = tempF;
                    tempF = @(x,y,u,ux,uy,time) self.coeffFromSymbolicEquation(jfD,1,self.dimension,x,y,u,ux,uy,time);
                    self.coeff{end}.jdFunction = tempF;                    
                end
                if isempty(self.coeff{end}.fiFunction) && ~isempty(mfF)
                    tempF = @(x,y,u,ux,uy,time) self.coeffFromSymbolicEquation(mfF,1,self.dimension,x,y,u,ux,uy,time);
                    self.coeff{end}.fiFunction = tempF;                  
                end
            end
            self.coeff{end}.symbolicFunction = ps.Results.symbolicFunction; 
            self.coeff{end}.parameters = ps.Results.parameters;
            self.coeff{end}.cFunction = [];
            self.coeff{end}.aFunction = [];
            self.coeff{end}.dFunction = [];
            self.coeff{end}.fFunction = [];
            if ~isempty(self.coeff{end}.symbolicFunction)                
                [cFunction,aFunction,fFunction,dFunction] = pdeFunctions('equationsFunction',self.coeff{end}.symbolicFunction,'forceGen',ps.Results.forceGen);
                self.coeff{end}.cFunction = @(p,t,u,time) cFunction(ps.Results.parameters,p,t,u,time);
                self.coeff{end}.aFunction = @(p,t,u,time) aFunction(ps.Results.parameters,p,t,u,time);
                self.coeff{end}.fFunction = @(p,t,u,time) fFunction(ps.Results.parameters,p,t,u,time);
                self.coeff{end}.dFunction = @(p,t,u,time) dFunction(ps.Results.parameters,p,t,u,time);                
            end
        end
        function ff = get.fFunction(self)
            ff = @self.fFunction_impl;
        end
        function cf = get.cFunction(self)
            cf = @self.cFunction_impl;
        end   
        function df = get.dFunction(self)
            df = @self.dFunction_impl;
        end   
        function af = get.aFunction(self)
            af = @self.aFunction_impl;
        end       
        function [cMatrix,aMatrix,dMatrix,F] = getMatrices(self,varargin)
            %getMatrices Get area integral matrices as described in PDE Toolbox function assema
            %   [cMatrix,aMatrix,dMatrix,F] = getMatrices(PARAM1,VALUE1,PARAM2,VALUE2,...) where PARAM can be:
            %       'u'   : Output variable of PDE as an array (#Dimensions X #Nodes) long.
            %       'time': Value of time as double.
            %       'p','e','t': Values of externally specified mesh of structure as described in initmesh.
            %
            %   Example: Get cMatrix,aMatrix,F matrices corresponding to u = u0 for coeffsObject, coeffs.
            %       [cMatrix,aMatrix,F] = coeffs.getMatrices('u',u0);
            %
            %   See also assema initmesh meshObject.
            ps = self.parser.parseInputArgs(varargin{:});    
            u = ps.Results.u;
            time = ps.Results.time;
            p = ps.Results.p;
            t = ps.Results.t;
            if isempty(p)
                p = self.geometry.mesh.p;
            end
            if isempty(t)
                t = self.geometry.mesh.t;
            end
            dMatrix = [];
            if isempty(u) && isempty(time)
                [cMatrix,aMatrix,F] = assema(p,t,self.cFunction,self.aFunction,self.fFunction);
                if ~isempty(self.dFunction)
                    [~,dMatrix] = assema(p,t,[],self.dFunction,[]);
                end
            elseif isempty(time)
                [cMatrix,aMatrix,F] = assema(p,t,self.cFunction,self.aFunction,self.fFunction,u);
                if ~isempty(self.dFunction)
                    [~,dMatrix] = assema(p,t,[],self.dFunction,[],u);
                end
            elseif isempty(u)
                [cMatrix,aMatrix,F] = assema(p,t,self.cFunction,self.aFunction,self.fFunction,time);
                if ~isempty(self.dFunction)
                    [~,dMatrix] = assema(p,t,[],self.dFunction,[],time);
                end
            else
                [cMatrix,aMatrix,F] = assema(p,t,self.cFunction,self.aFunction,self.fFunction,u,time);
                if ~isempty(self.dFunction)
                    [~,dMatrix] = assema(p,t,[],self.dFunction,[],u,time);
                end
            end
        end
        function [cMatrix,aMatrix,dMatrix,F] = getJacobian(self,varargin)
           %getJacobian Get Jacobian of area integral matrices
            %   [K,M,F] = getJacobian(PARAM1,VALUE1,PARAM2,VALUE2,...) where PARAM can be:
            %       'u'   : Output variable of PDE as an array (#Dimensions X #Nodes) long.
            %       'time': Value of time as double.
            %       'p','e','t': Values of externally specified mesh of structure as described in initmesh.
            %   getJacobian is applicable only for the case when the 'symbolicEquationFunction'
            %   input is used and is useful only when both the PDE is nonlinear and a custom
            %   nonlinear solver is used.
            %        
            ps = self.parser.parseInputArgs(varargin{:});
            u = ps.Results.u;
            time = ps.Results.time;
            p = ps.Results.p;
            t = ps.Results.t;
            if isempty(p)
                p = self.geometry.mesh.p;
            end
            if isempty(t)
                t = self.geometry.mesh.t;
            end
            dMatrix = [];
            if isempty(u) && isempty(time)
                [cMatrix,aMatrix,F] = assema(p,t,@self.jcFunction,@self.jaFunction,@self.fFunction);
            elseif isempty(time)
                [cMatrix,aMatrix,F] = assema(p,t,@self.jcFunction,@self.jaFunction,@self.fFunction,u);
            elseif isempty(u)
                [cMatrix,aMatrix,F] = assema(p,t,@self.jcFunction,@self.jaFunction,@self.fFunction,time);
            else
                [cMatrix,aMatrix,F] = assema(p,t,@self.jcFunction,@self.jaFunction,@self.fFunction,u,time);
            end
        end
    end
    methods(Hidden = true)
        function [c] = jcFunction(self,p,t,u,time)
            c = self.cGeneralFunction_impl(p,t,u,time,false);
        end
        function [a] = jaFunction(self,p,t,u,time)
            a = self.aGeneralFunction_impl(p,t,u,time,false);
        end          
    end
    methods(Access = private, Hidden = true)
        function [f] = fFunction_impl(self,p,t,u,time)
            regionToIds = self.geometry.regionToId;
            N = self.dimension;
            nt = size(t,2);
            f = zeros(N,nt);
            % Triangle point indices
            it1=t(self.geometry.mesh.triangleFirstPointIndex,:);
            it2=t(self.geometry.mesh.triangleSecondPointIndex,:);
            it3=t(self.geometry.mesh.triangleThirdPointIndex,:);
            % Find centroids of triangles
            xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
            ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;
            if ~isempty(u)
                uintrp = pdeintrp(p,t,u); % Interpolated values at centroids
                uintrp = reshape(uintrp,[],N); % matrix with N column
                uintrp = uintrp'; % change to row vectors
                [ux,uy] = pdegrad(p,t,u); % Approximate derivatives
            else
                uintrp = [];
                ux = [];
                uy = [];
            end
            for cix=1:length(self.coeff)
                if ~isempty(self.coeff{cix}.fFunction)
                    regionId = regionToIds(self.coeff{cix}.region);
                    ix = find(t(self.geometry.mesh.triangleRegionIndex,:) == regionId);
                    f(:,ix) = self.coeff{cix}.fFunction(p,t(:,ix),u,time);
                    continue;
                end
                fiFunction = self.coeff{cix}.fiFunction;
                if ~isempty(fiFunction)
                    regionId = regionToIds(self.coeff{cix}.region);
                    ix = find(t(self.geometry.mesh.triangleRegionIndex,:) == regionId);
                    if ix
                        xc = xpts(ix);
                        yc = ypts(ix);
                        if ~isempty(u)
                            uc = uintrp(:,ix);
                            uxc = ux(:,ix);
                            uyc = uy(:,ix);
                        end
                        for k=1:length(xc)
                            if ~isempty(u)
                                if ~self.triangleIndividualCoordinatesForm
                                    fi = fiFunction(xc(k),yc(k),uc(:,k),uxc(:,k),uyc(:,k),time);
                                else
                                    p1 = p(:,t(self.geometry.mesh.triangleFirstPointIndex,ix(k)));
                                    p2 = p(:,t(self.geometry.mesh.triangleSecondPointIndex,ix(k)));
                                    p3 = p(:,t(self.geometry.mesh.triangleThirdPointIndex,ix(k)));
                                    fi = fiFunction([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)],uc(:,k),uxc(:,k),uyc(:,k),time);
                                end
                            else
                                if ~self.triangleIndividualCoordinatesForm
                                    fi = fiFunction(xc(k),yc(k),[],[],[],time);
                                else
                                    p1 = p(:,t(self.geometry.mesh.triangleFirstPointIndex,ix(k)));
                                    p2 = p(:,t(self.geometry.mesh.triangleSecondPointIndex,ix(k)));
                                    p3 = p(:,t(self.geometry.mesh.triangleThirdPointIndex,ix(k)));
                                    fi = fiFunction([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)],[],[],[],time);
                                end
                            end
                            for n=1:N
                                f(n,ix(k)) = fi(n);
                            end
                        end
                    end
                end
            end
        end
        function [c] = cFunction_impl(self,p,t,u,time)
            c = self.cGeneralFunction_impl(p,t,u,time,true);
        end
        function [d] = dFunction_impl(self,p,t,u,time)
            regionToIds = self.geometry.regionToId;
            N = self.dimension;
            nt = size(t,2);
            d = zeros(N*N,nt);
            % Triangle point indices
            it1=t(self.geometry.mesh.triangleFirstPointIndex,:);
            it2=t(self.geometry.mesh.triangleSecondPointIndex,:);
            it3=t(self.geometry.mesh.triangleThirdPointIndex,:);
            % Find centroids of triangles
            xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
            ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;
            if ~isempty(u)
                uintrp = pdeintrp(p,t,u); % Interpolated values at centroids
                uintrp = reshape(uintrp,[],N); % matrix with N column
                uintrp = uintrp'; % change to row vectors
                [ux,uy] = pdegrad(p,t,u); % Approximate derivatives
            else
                uintrp = [];
                ux = [];
                uy = [];
            end
            for cix=1:length(self.coeff)
                if ~isempty(self.coeff{cix}.dFunction)
                    regionId = regionToIds(self.coeff{cix}.region);
                    ix = find(t(self.geometry.mesh.triangleRegionIndex,:) == regionId);
                    d(:,ix) = self.coeff{cix}.dFunction(p,t(:,ix),u,time);
                    continue;
                end
                dijFunction = self.coeff{cix}.dijFunction;
                if ~isempty(dijFunction)
                    regionId = regionToIds(self.coeff{cix}.region);
                    ix = find(t(self.geometry.mesh.triangleRegionIndex,:) == regionId);
                    if ix
                        xc = xpts(ix);
                        yc = ypts(ix);
                        if ~isempty(u)
                            uc = uintrp(:,ix);
                            uxc = ux(:,ix);
                            uyc = uy(:,ix);
                        end
                        for k=1:length(xc)
                            if ~isempty(u)
                                if ~self.triangleIndividualCoordinatesForm
                                    dij = dijFunction(xc(k),yc(k),uc(:,k),uxc(:,k),uyc(:,k),time);
                                else
                                    p1 = p(:,t(self.geometry.mesh.triangleFirstPointIndex,ix(k)));
                                    p2 = p(:,t(self.geometry.mesh.triangleSecondPointIndex,ix(k)));
                                    p3 = p(:,t(self.geometry.mesh.triangleThirdPointIndex,ix(k)));
                                    dij = dijFunction([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)],uc(:,k),uxc(:,k),uyc(:,k),time);
                                end
                            else
                                if ~self.triangleIndividualCoordinatesForm
                                    dij = dijFunction(xc(k),yc(k),[],[],[],time);
                                else
                                    p1 = p(:,t(self.geometry.mesh.triangleFirstPointIndex,ix(k)));
                                    p2 = p(:,t(self.geometry.mesh.triangleSecondPointIndex,ix(k)));
                                    p3 = p(:,t(self.geometry.mesh.triangleThirdPointIndex,ix(k)));
                                    dij = dijFunction([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)],[],[],[],time);
                                end
                            end
                            d(:,ix(k)) = dij(:);
                        end
                    end
                end
            end
        end
        function [c] = cGeneralFunction_impl(self,p,t,u,time,iscFunction)
            regionToIds = self.geometry.regionToId;
            N = self.dimension;
            nt = size(t,2);
            c = zeros(4*N^2,nt);
            % Triangle point indices
            it1=t(self.geometry.mesh.triangleFirstPointIndex,:);
            it2=t(self.geometry.mesh.triangleSecondPointIndex,:);
            it3=t(self.geometry.mesh.triangleThirdPointIndex,:);
            % Find centroids of triangles
            xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
            ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;
            if ~isempty(u)
                uintrp = pdeintrp(p,t,u); % Interpolated values at centroids
                uintrp = reshape(uintrp,[],N); % matrix with N column
                uintrp = uintrp'; % change to row vectors
                [ux,uy] = pdegrad(p,t,u); % Approximate derivatives
            else
                uintrp = [];
                ux = [];
                uy = [];
            end            
            for cix=1:length(self.coeff)
                if ~isempty(self.coeff{cix}.cFunction)
                    regionId = regionToIds(self.coeff{cix}.region);
                    ix = find(t(self.geometry.mesh.triangleRegionIndex,:) == regionId);
                    c(:,ix) = self.coeff{cix}.cFunction(p,t(:,ix),u,time);
                    continue;
                end
                cijFunction = self.coeff{cix}.cijFunction;
                if ~iscFunction
                    cijFunction = self.coeff{cix}.jcFunction;
                end
                if ~isempty(cijFunction)
                    regionId = regionToIds(self.coeff{cix}.region);
                    ix = find(t(self.geometry.mesh.triangleRegionIndex,:) == regionId);
                    if ix
                        for k=ix
                            if ~isempty(u)
                                cij = cijFunction(xpts(k),ypts(k),uintrp(:,k),ux(:,k),uy(:,k),time);
                            else
                                cij = cijFunction(xpts(k),ypts(k),[],[],[],time);
                            end
                            count = 1;
                            for n=1:N
                                for m=1:N
                                    c(count,k) = cij(2*(m-1)+1, 2*(n-1)+1);
                                    count = count + 1;
                                    c(count,k) = cij(2*(m-1)+2, 2*(n-1)+1);
                                    count = count + 1;
                                    c(count,k) = cij(2*(m-1)+1, 2*(n-1)+2);
                                    count = count + 1;
                                    c(count,k) = cij(2*(m-1)+2, 2*(n-1)+2);
                                    count = count+1;
                                end
                            end
                        end
                    end
                end
            end
        end
        function [a] = aFunction_impl(self,p,t,u,time)
            a = self.aGeneralFunction_impl(p,t,u,time,true);
        end
        function [a] = aGeneralFunction_impl(self,p,t,u,time,isaFunction)
            regionToIds = self.geometry.regionToId;
            N = self.dimension;
            nt = size(t,2);
            a = zeros(N*N,nt);
            % Triangle point indices
            it1=t(self.geometry.mesh.triangleFirstPointIndex,:);
            it2=t(self.geometry.mesh.triangleSecondPointIndex,:);
            it3=t(self.geometry.mesh.triangleThirdPointIndex,:);
            % Find centroids of triangles
            xpts=(p(1,it1)+p(1,it2)+p(1,it3))/3;
            ypts=(p(2,it1)+p(2,it2)+p(2,it3))/3;
            if ~isempty(u)
                uintrp = pdeintrp(p,t,u); % Interpolated values at centroids
                uintrp = reshape(uintrp,[],N); % matrix with N column
                uintrp = uintrp'; % change to row vectors
                [ux,uy] = pdegrad(p,t,u); % Approximate derivatives
            else
                uintrp = [];
                ux = [];
                uy = [];
            end
            for cix=1:length(self.coeff)
                if ~isempty(self.coeff{cix}.aFunction)
                    regionId = regionToIds(self.coeff{cix}.region);
                    ix = find(t(self.geometry.mesh.triangleRegionIndex,:) == regionId);
                    a(:,ix) = self.coeff{cix}.aFunction(p,t(:,ix),u,time);
                    continue;
                end
                aijFunction = self.coeff{cix}.aijFunction;
                if ~isaFunction
                    aijFunction = self.coeff{cix}.jaFunction;
                end
                if ~isempty(aijFunction)
                    regionId = regionToIds(self.coeff{cix}.region);
                    ix = find(t(self.geometry.mesh.triangleRegionIndex,:) == regionId);
                    if ix
                        xc = xpts(ix);
                        yc = ypts(ix);
                        if ~isempty(u)
                            uc = uintrp(:,ix);
                            uxc = ux(:,ix);
                            uyc = uy(:,ix);
                        end
                        for k=1:length(xc)
                            if ~isempty(u)
                                if ~self.triangleIndividualCoordinatesForm
                                    aij = aijFunction(xc(k),yc(k),uc(:,k),uxc(:,k),uyc(:,k),time);
                                else
                                    p1 = p(:,t(self.geometry.mesh.triangleFirstPointIndex,ix(k)));
                                    p2 = p(:,t(self.geometry.mesh.triangleSecondPointIndex,ix(k)));
                                    p3 = p(:,t(self.geometry.mesh.triangleThirdPointIndex,ix(k)));
                                    aij = aijFunction([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)],uc(:,k),uxc(:,k),uyc(:,k),time);
                                end
                            else
                                if ~self.triangleIndividualCoordinatesForm
                                    aij = aijFunction(xc(k),yc(k),[],[],[],time);
                                else
                                    p1 = p(:,t(self.geometry.mesh.triangleFirstPointIndex,ix(k)));
                                    p2 = p(:,t(self.geometry.mesh.triangleSecondPointIndex,ix(k)));
                                    p3 = p(:,t(self.geometry.mesh.triangleThirdPointIndex,ix(k)));
                                    aij = aijFunction([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)],[],[],[],time);
                                end
                            end
                            a(:,ix(k)) = aij(:);
                        end
                    end
                end
            end
        end
    end
    methods(Static, Hidden = true)        
        function [cij] = cijConstantDiagFunctionHelper(cVal,N)
            cij = sparse(2*N,2*N);
            if ~isscalar(cVal)
                for n=1:size(cij,1)
                    cij(n,n) = cVal(n);
                end
            else
                for n=1:size(cij,1)
                    cij(n,n) = cVal;
                end
            end
        end
        function fi = fiConstantFunctionHelper(fVal,N)
            if isscalar(fVal)
                fi = fVal*ones(N,1);
            else
                fi = fVal;
            end
        end
        function [aij] = aijConstantDiagFunctionHelper(aVal,N)
            aij = sparse(N,N);
            if ~isscalar(aVal)
                for n=1:size(aij,1)
                    aij(n,n) = aVal(n);
                end
            else
                for n=1:size(aij,1)
                    aij(n,n) = aVal;
                end
            end
        end
    end
    methods(Static,Access = private)
        function fi = zeroFFunction(N)
            fi = zeros(N,1);
        end
        function cij = coeffFromSymbolicEquation(mfC,intFactor,N,x,y,u,ux,uy,time)
            if isempty(time)
                time = 0;
            end
            uxy = zeros(2*N,1);
            if isempty(u)
                u = zeros(N,1);
            else
                uxy(1:2:end-1,1) = ux;
                uxy(2:2:end,1) = uy;
            end
            ucell = num2cell(u);
            uxycell = num2cell(uxy);
            cij = mfC(x,y,time,ucell{:},uxycell{:});
%             switch length(u)
%                 case 1                     
%                     if isempty(u)
%                         cij = mfC(x,y,time,0,0,0);
%                     else
%                         cij = mfC(x,y,time,u(1),ux(1),uy(1));
%                     end
%                 case 2
%                     if isempty(u)
%                         cij = mfC(x,y,time,0,0,0,0,0,0);
%                     else
%                         cij = mfC(x,y,time,u(1),u(2),ux(1),uy(1),ux(2),uy(2));
%                     end
%                 case 3                    
%                     if isempty(u)
%                         cij = mfC(x,y,time,0,0,0,0,0,0,0,0,0);
%                     else
%                         cij = mfC(x,y,time,u(1),u(2),u(3),ux(1),uy(1),ux(2),uy(2),ux(3),uy(3));
%                     end
%                 case 4                    
%                     if isempty(u)
%                         cij = mfC(x,y,time,0,0,0,0,0,0,0,0,0,0,0,0);
%                     else
%                         cij = mfC(x,y,time,u(1),u(2),u(3),u(4),ux(1),uy(1),ux(2),uy(2),ux(3),uy(3),ux(4),uy(4));
%                     end                    
%                 case 5                    
%                     if isempty(u)
%                         cij = mfC(x,y,time,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
%                     else
%                         cij = mfC(x,y,time,u(1),u(2),u(3),u(4),u(5),ux(1),uy(1),ux(2),uy(2),ux(3),uy(3),ux(4),uy(4),ux(5),uy(5));
%                     end
%             end
            cij = reshape(cij,intFactor*N,[]);
        end
        function [mfC,mfA,mfF,mfD,jfC,jfA,jfD] = generateMatrixFunctionsFromSymbolic(varargin)
            equations = [];
            variables = [];
            index = 1;
            displayCoefficients = true;
            for k=1:2:length(varargin)
                a = varargin(k);
                b = varargin(k+1);
                if strcmp(a,'equations')
                    equations = b{1};
                elseif strcmp(a,'variables')
                    variables = b{1};
                elseif strcmp(a,'index')
                    index = b{1};                     
                elseif strcmp(a,'displayCoefficients')
                    displayCoefficients = b{1}; 
                end
            end
            uprime = sym(zeros(1,length(variables)*2));
            ddt = sym(zeros(1,length(variables)));
            d2d2t = sym(zeros(1,length(variables)));
            for k=1:length(variables)
                uprime(2*k-1) = sym(strcat('d',char(variables(k)),'dx'));
                uprime(2*k) = sym(strcat('d',char(variables(k)),'dy'));
                ddt(k) = sym(strcat('d',char(variables(k)),'dt'));
                d2d2t(k) = sym(strcat('d2',char(variables(k)),'d2t'));
            end
            variablesPlus = [sym('x') sym('y') sym('time') variables uprime];
            % extract ddt* and d2dt* terms
            nond2d2tequations = subs(equations,d2d2t,0*d2d2t);
            Dequations = equations - nond2d2tequations;
            equations = nond2d2tequations;
            if ~isempty(find(Dequations ~= 0,1))
                D = pdetbplus.coeffsObject.extractDMatrix(Dequations,d2d2t,displayCoefficients);
                filename = strcat('Dfunction',num2str(index));
                mfD = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(D,variablesPlus,filename);
                jD = jacobian(D*d2d2t.',d2d2t);
                filename = strcat('jDfunction',num2str(index));
                jfD = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(jD,variablesPlus,filename);
            else
                nonddtequations = subs(equations,ddt,0*ddt);
                Dequations = equations - nonddtequations;
                equations = nonddtequations;                
                if ~isempty(find(Dequations ~= 0,1))
                    D = pdetbplus.coeffsObject.extractDMatrix(Dequations,ddt,displayCoefficients);
                    filename = strcat('Dfunction',num2str(index));
                    mfD = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(D,variablesPlus,filename);
                    jD = jacobian(D*ddt.',ddt);
                    filename = strcat('jDfunction',num2str(index));
                    jfD = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(jD,variablesPlus,filename);
                else                    
                    mfD = [];
                    jfD = [];
                end
            end
            % extract ddx*, ddy* terms, C matrix
            nonCequations = subs(equations,[sym('ddx'),sym('ddy')],[0,0]);
            Cequations = -(equations - nonCequations); % "-" because the Div template expects this form
            C = pdetbplus.coeffsObject.extractCMatrix(Cequations,uprime,displayCoefficients);            
            filename = strcat('Cfunction',num2str(index));
            mfC = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(C,variablesPlus,filename);
            jC = jacobian(C*uprime.',uprime);
            filename = strcat('jCfunction',num2str(index));
            jfC = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(jC,variablesPlus,filename);
            % extract F vector
            Fequations = -subs(nonCequations,variables,0*variables); % "-" because "f" needs to be on the RHS            
            if displayCoefficients
                display Fvector; pretty(sym(Fequations));
            end
            if ~isempty(find(Fequations ~= 0,1))
                filename = strcat('Ffunction',num2str(index));
                mfF = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(Fequations,variablesPlus,filename);
            else
                mfF = [];
            end
            % extract A matrix
            Aequations = (nonCequations - Fequations);            
            if ~isempty(find(Aequations ~= 0,1))
                A = pdetbplus.coeffsObject.extractAMatrix(Aequations,variables,displayCoefficients);
                filename = strcat('Afunction',num2str(index));
                mfA = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(A,variablesPlus,filename);
                jA = jacobian(A*variables.',variables);
                filename = strcat('jAfunction',num2str(index));
                jfA = pdetbplus.coeffsObject.convertSymbolicMatrixToFunction(jA,variablesPlus,filename);
            else
                mfA = [];
                jfA = [];
            end
        end     
        function C = extractCMatrix(equations,uprime,displayCoefficients)
            M = equationsToMatrix(equations,[sym('ddx') sym('ddy')]);
            % M is a N x 2 matrix; we convert this to a 2N vector
            Mt = M.';
            Mvec = Mt(:);
            [Num,Den] = numden(Mvec);
            C = sym(zeros(length(Num),length(Num)));            
            for k=1:size(C,1) % loop over rows of M2
                for j=1:size(C,2) % loop over columns of M2
                    % get coeffs w.r.t. uprime; D contains the polynomial terms [1,x,x^2,...]; B contains
                    % the corresponding coefficients [a0,a1,a2,...]
                    [B,D] = coeffs(Num(k),uprime(j));
                    if ~isempty(B)
                        % set all components of D to zero where uprime(j) occurs. This
                        % is to exclude the term that does not contain uprime(j)
                        Dz = subs(D,uprime(j),0);
                        izx = find(Dz == 0);
                        if ~isempty(izx)
                            % assemble and sum the terms sans the term that does not
                            % contain uprime(j)
                            C(k,j) = simplify(D(izx)*B(izx)'/(uprime(j)*Den(k)));
                            % we've already accounted for uprime(j) in row k so set
                            % all terms that have uprime(j) to 0
                            Num(k) = subs(Num(k),uprime(j),0);
                        end
                    end
                end
            end
            if displayCoefficients
                display Cmatrix;
                pretty(C);
            end
        end     
        function A = extractAMatrix(equations,u,displayCoefficients)
            A = equationsToMatrix(equations,u);
            if displayCoefficients
                display Amatrix;
                pretty(A);
            end
        end
        function D = extractDMatrix(equations,u,displayCoefficients)
            D = equationsToMatrix(equations,u);
            if displayCoefficients
                display Dmatrix;
                pretty(D);
            end
        end
        function XFunction = convertSymbolicMatrixToFunction(X,variables,filename)
            XFunction = matlabFunction(X(:),'vars',variables,'file',filename);
        end
    end
end