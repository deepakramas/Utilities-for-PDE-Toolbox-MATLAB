classdef xyFunctionClass
    properties
        sol;
        mesh;
        element;
    end
    methods
        function self = xyFunctionClass(sol,mesh)
            self.sol = sol;
            self.mesh = mesh;
            self.element.minx = nan*ones(size(self.mesh.t,2),1);
            self.element.miny = nan*ones(size(self.mesh.t,2),1);
            self.element.maxx = nan*ones(size(self.mesh.t,2),1);
            self.element.maxy = nan*ones(size(self.mesh.t,2),1);
            self.element.pt1Coeff = cell(size(self.mesh.t,2),1);
            self.element.pt2Coeff = cell(size(self.mesh.t,2),1);
            self.element.pt3Coeff = cell(size(self.mesh.t,2),1);
            for k=1:length(self.mesh.t)
                a1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleFirstPointIndex,k));
                a2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleFirstPointIndex,k));
                b1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleSecondPointIndex,k));
                b2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleSecondPointIndex,k));
                c1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleThirdPointIndex,k));
                c2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleThirdPointIndex,k));
                self.element.minx(k) = min([a1 b1 c1]);
                self.element.miny(k) = min([a2 b2 c2]);
                self.element.maxx(k) = max([a1 b1 c1]);
                self.element.maxy(k) = max([a2 b2 c2]);
                T1 = a2*c1 - a2*b1 + b1*c2 - b2*c1 + a1*(b2 - c2);
                self.element.pt1Coeff{k} =  1/T1*[(b2 - c2); -(b1 - c1); (b1*c2 - b2*c1)];
                self.element.pt2Coeff{k} = 1/T1*[-(a2 - c2); (a1 - c1); -(a1*c2 - a2*c1)];
                self.element.pt3Coeff{k} = 1/T1*[(a2 - b2); -(a1 - b1); (a1*b2 - a2*b1)];
            end
        end
        function xySol = xyFunction(self,xvec,yvec)
            if isscalar(yvec)&& ~isscalar(xvec)
                yvec = yvec*ones(size(xvec));
            elseif isscalar(xvec)&& ~isscalar(yvec)
                xvec = xvec*ones(size(yvec));
            end
            if isrow(xvec)
                xvec = xvec.';
            end
            xySol = nan*ones(size(xvec,1),size(self.sol,2));  % always pts along row and other dimension is for other solution values
            for ix=1:size(xvec,1)
                x = xvec(ix);
                y = yvec(ix);
                candidateElementIndices = find(x >= self.element.minx & x <= self.element.maxx & y >= self.element.miny & y <= self.element.maxy);
                for kix = candidateElementIndices'
                    a1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleFirstPointIndex,kix));
                    a2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleFirstPointIndex,kix));
                    b1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleSecondPointIndex,kix));
                    b2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleSecondPointIndex,kix));
                    c1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleThirdPointIndex,kix));
                    c2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleThirdPointIndex,kix));
                    % Based on reference: http://blogs.msdn.com/b/rezanour/archive/2011/08/07/barycentric-coordinates-and-point-in-triangle-tests.aspx
                    % barycentric variables
                    u = [b1 b2 0]' - [a1 a2 0]';
                    v = [c1 c2 0]' - [a1 a2 0]';
                    w = [x y 0]' - [a1 a2 0]';
                    % cross products that define the barycentric
                    % coordinates,r and t
                    vXw = cross(v, w); vXu = cross(v, u);
                    uXw = cross(u, w); uXv = cross(u, v);
                    r = vXw(3)/vXu(3);
                    t = uXw(3)/uXv(3);
                    inside = (r >= 0 && r <= 1) && (t >= 0 && t <= 1) && (r + t) <= 1;
                    if inside % we've found our triangle
                        sola = self.sol(self.mesh.t(self.mesh.triangleFirstPointIndex,kix),:);
                        solb = self.sol(self.mesh.t(self.mesh.triangleSecondPointIndex,kix),:);
                        solc = self.sol(self.mesh.t(self.mesh.triangleThirdPointIndex,kix),:);
                        xySol(ix,:) = [x y 1]*(self.element.pt1Coeff{kix}*sola + self.element.pt2Coeff{kix}*solb + self.element.pt3Coeff{kix}*solc);
                        break;
                    end
                end
            end
            % transform it back to row if input was row
            if isrow(yvec)
                xySol = xySol.';
            end
        end
        function xyToNode = elementIndexFunction(self,xvec,yvec)
            if isscalar(yvec)&& ~isscalar(xvec)
                yvec = yvec*ones(size(xvec));
            elseif isscalar(xvec)&& ~isscalar(yvec)
                xvec = xvec*ones(size(yvec));
            end
            if isrow(xvec)
                xvec = xvec.';
            end
            xyToNode = nan*ones(size(xvec,1),6);
            for ix=1:size(xvec,1)
                x = xvec(ix);
                y = yvec(ix);
                candidateElementIndices = find(x >= self.element.minx & x <= self.element.maxx & y >= self.element.miny & y <= self.element.maxy);
                for kix = candidateElementIndices'
                    a1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleFirstPointIndex,kix));
                    a2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleFirstPointIndex,kix));
                    b1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleSecondPointIndex,kix));
                    b2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleSecondPointIndex,kix));
                    c1 = self.mesh.p(1,self.mesh.t(self.mesh.triangleThirdPointIndex,kix));
                    c2 = self.mesh.p(2,self.mesh.t(self.mesh.triangleThirdPointIndex,kix));
                    % Based on reference: http://blogs.msdn.com/b/rezanour/archive/2011/08/07/barycentric-coordinates-and-point-in-triangle-tests.aspx
                    % barycentric variables
                    u = [b1 b2 0]' - [a1 a2 0]';
                    v = [c1 c2 0]' - [a1 a2 0]';
                    w = [x y 0]' - [a1 a2 0]';
                    % cross products that define the barycentric
                    % coordinates,r and t
                    vXw = cross(v, w); vXu = cross(v, u);
                    uXw = cross(u, w); uXv = cross(u, v);
                    r = vXw(3)/vXu(3);
                    t = uXw(3)/uXv(3);
                    inside = (r >= (-100*eps) && r <= (1+100*eps)) && (t >= (-100*eps) && t <= (1+100*eps)) && (r + t) <= (1+100*eps);
                    if inside % we've found our triangle
                        node1 = self.mesh.t(self.mesh.triangleFirstPointIndex,kix);
                        node2 = self.mesh.t(self.mesh.triangleSecondPointIndex,kix);
                        node3 = self.mesh.t(self.mesh.triangleThirdPointIndex,kix);
                        nodalContributions = [x y 1]*[self.element.pt1Coeff{kix}, self.element.pt2Coeff{kix}, self.element.pt3Coeff{kix}];
                        xyToNode(ix,:) = [node1,node2,node3,nodalContributions];
                        break;
                    end
                end
            end
            % transform it back to row if input was row
            if isrow(yvec)
                xyToNode = xyToNode.';
            end
        end
        function [dudxContrib,dudyContrib] = gradientContributions(self)
            dudxContrib = sparse(size(self.mesh.t,2),size(self.mesh.p,2));
            dudyContrib = sparse(size(self.mesh.t,2),size(self.mesh.p,2));
            for kix = 1:size(self.mesh.t,2)
                dudxContrib(kix,self.mesh.t(1:3,kix)) = [self.element.pt1Coeff{kix}(1), self.element.pt2Coeff{kix}(1), self.element.pt3Coeff{kix}(1)];
                dudyContrib(kix,self.mesh.t(1:3,kix)) = [self.element.pt1Coeff{kix}(2), self.element.pt2Coeff{kix}(2), self.element.pt3Coeff{kix}(2)];
            end
        end
    end
end