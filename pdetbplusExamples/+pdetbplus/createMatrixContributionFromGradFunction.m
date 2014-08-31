function P = createMatrixContributionFromGradFunction(p,t,fFunction1,regionIndex)
% create node to triangles map
nodeToTriangles = cell(size(p,2),1);
if nargin < 5
    regionSpecified = false;
else
    regionSpecified = true;
end
for k=1:size(t,2)
    if (regionSpecified && t(5,k) == regionIndex) || ~regionSpecified
        nodeToTriangles{t(1,k)} = [nodeToTriangles{t(1,k)} k];
        nodeToTriangles{t(2,k)} = [nodeToTriangles{t(2,k)} k];
        nodeToTriangles{t(3,k)} = [nodeToTriangles{t(3,k)} k];
    end
end
% triangle coordinates : 1 (a1,a2) 2 (b1,b2) 3 (c1,c2)
% ux = (b2 - c2)/(a2*c1 - a2*b1 + b1*c2 - b2*c1 + a1*(b2 - c2)) * u1 +
% -(a2 - c2)/(a2*c1 - a2*b1 + b1*c2 - b2*c1 + a1*(b2 - c2))*u2 +
% (a2 - b2)/(a2*c1 - a2*b1 + b1*c2 - b2*c1 + a1*(b2 - c2))*u3

% uy = (b1 - c1)/(a1*c2 - a1*b2 - b1*c2 + b2*c1 + a2*(b1 - c1))*u1 +
% (a1 - c1)/(a1*b2 - a2*b1 - a1*c2 + a2*c1 + b1*c2 - b2*c1)*u2 +
% -(a1 - b1)/(a1*b2 - a2*b1 - a1*c2 + a2*c1 + b1*c2 - b2*c1)*u3

N = 1;
    function [f] = fFunction2(p,t,~,~)       
        nt = size(t,2);
        f = zeros(N,nt);
        triangles = nodeToTriangles{nodeIndex};
        for tix = triangles
            node1Index = t(1,tix);
            a1 = p(1,node1Index);
            a2 = p(2,node1Index);
            node2Index = t(2,tix);
            b1 = p(1,node2Index);
            b2 = p(2,node2Index);
            node3Index = t(3,tix);
            c1 = p(1,node3Index);
            c2 = p(2,node3Index);
            [uxcoeff,uycoeff] = fFunction1((a1+b1+c1)/3,(a2+b2+c2)/3);
            switch nodeIndex
                case node1Index                    
                    f(N,tix) = uxcoeff*((b2 - c2)/(a2*c1 - a2*b1 + b1*c2 - b2*c1 + a1*(b2 - c2))) + ...
                        uycoeff*((b1 - c1)/(a1*c2 - a1*b2 - b1*c2 + b2*c1 + a2*(b1 - c1)));
                case node2Index
                    f(N,tix) = uxcoeff*(-(a2 - c2)/(a2*c1 - a2*b1 + b1*c2 - b2*c1 + a1*(b2 - c2))) + ...
                        uycoeff*((a1 - c1)/(a1*b2 - a2*b1 - a1*c2 + a2*c1 + b1*c2 - b2*c1));
                case node3Index
                    f(N,tix) = uxcoeff*((a2 - b2)/(a2*c1 - a2*b1 + b1*c2 - b2*c1 + a1*(b2 - c2))) + ...
                        uycoeff*(-(a1 - b1)/(a1*b2 - a2*b1 - a1*c2 + a2*c1 + b1*c2 - b2*c1));
            end           
        end       
    end
P = zeros(N*size(p,2),N*size(p,2));
for nodeIndex=1:size(p,2)
    [~,~,F] = assema(p,t,0,0,@fFunction2);
    P(:,nodeIndex) = F;
end
P = sparse(P);
end