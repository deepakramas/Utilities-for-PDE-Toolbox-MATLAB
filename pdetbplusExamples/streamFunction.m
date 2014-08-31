function streamFunction
clf;hold on;
import pdetbplus.*;

op = {pointObject(-1,-1),pointObject(-1,0),pointObject(0,0),pointObject(0,1),pointObject(1,1),pointObject(1,-1)};
square = geometryObject.createPolygon('name','square','points',op,'leftRegion','cavity','rightRegion','material');
% ip = {pointObject(0.3,-0.5),pointObject(0.3,-0.3),pointObject(0.7,-0.3),pointObject(0.7,-0.5)};
% hole = geometryObject.createPolygon('name','hole','points',ip,'leftRegion','material','rightRegion','cavity');
ip = {pointObject(0.7,-0.5),pointObject(0.7,-0.3),pointObject(0.3,-0.3),pointObject(0.3,-0.5)};
hole = geometryObject.createPolygon('name','hole','points',ip,'leftRegion','cavity','rightRegion','material');

squareWithHole = square + hole;
squareWithHole.exteriorRegion = 'cavity';
% squareWithHole = squareWithHole.rotate(pointObject(0,0),pi/6);
squareWithHole = squareWithHole.initMesh('showMesh',false,'Hmax',1/20);

coeff = coeffsObject(squareWithHole,1);
cValue = 1;
coeff.add('region','material','aConstantValue',0,'cConstantValue',cValue,'fConstantValue',10);

bc = boundaryConditionObject(squareWithHole,1);
bc.add('name',square.boundary{4}.name,'dirichlet',1,'dimension',1);
bc.add('name',square.boundary{6}.name,'dirichlet',1000000,'dimension',1);

[Q,G,H,R] = bc.getMatrices();
[K,M,F] = coeff.getMatrices();
u = assempde(K,M,F,Q,G,H,R);
[cdudx,cdudy] = pdegrad(squareWithHole.mesh.p,squareWithHole.mesh.t,u);
[boundaryEdges,boundaryEdgeElements] = squareWithHole.getBoundaryEdges('cavity','material');
edgeToElement = zeros(max(boundaryEdges),1);
edgeToElement(boundaryEdges) = boundaryEdgeElements;

bc2 = boundaryConditionObject(squareWithHole,1);
    function [hval,rval,qval,gval] = bcPotentialToStream(x,y,u,t,edgeIndex)
        rval = []; hval = [];
        qval = 0;
        p1 = squareWithHole.mesh.e(1,edgeIndex);
        p2 = squareWithHole.mesh.e(2,edgeIndex);
        xy1 = squareWithHole.mesh.p(1:2,p1);
        xy2 = squareWithHole.mesh.p(1:2,p2);
        dp = xy2 - xy1;
        lengthEdge = norm(dp);
        nx = -dp(2)/lengthEdge;
        ny = dp(1)/lengthEdge;
        if (squareWithHole.mesh.e(6,edgeIndex) == squareWithHole.getRegionToId('material')) && (squareWithHole.mesh.e(7,edgeIndex) == squareWithHole.getRegionToId('cavity'))
            % nx, ny points from outer region to inner region. we need to flip the direction
            nx = -nx;
            ny = -ny;
        end
        cdudxVal = cdudx(edgeToElement(edgeIndex));
        cdudyVal = cdudy(edgeToElement(edgeIndex));
        gval = (-nx*cdudyVal + ny*cdudxVal);
    end

bc2.add('innerRegion','material','outerRegion','cavity','xyutFunction',@bcPotentialToStream);
[Q,G,H,R] = bc2.getMatrices();
boundaryNodes = squareWithHole.getBoundaryNodes('region1','material','region2','cavity');
[Hadd,Radd] = bc.getGlobalHR('nodes',boundaryNodes(2),'type','point','value',0,'dimension',1);
H = [H;Hadd];
R = [R;Radd];
v = assempde(K,M,F,Q,G,H,R);
pdeplot(squareWithHole.mesh.p,squareWithHole.mesh.e,squareWithHole.mesh.t,'xydata',u,'contour','on','mesh','off','levels',30);hold on;
pdeplot(squareWithHole.mesh.p,squareWithHole.mesh.e,squareWithHole.mesh.t,'xydata',v,'contour','on','mesh','off','levels',30);
shg;
end