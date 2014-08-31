function [u,converged] = solveGeomNonlinear(varargin)
uinitialPDE = [];
uoPDE = [];
coeff = [];
bc = [];
UpdateGeometry = false;
Hadd = [];
Radd = [];
Qadd = [];
Gadd = [];
ResidualAbsTol = 1e-3;
SolRelTol = 1e-3;
MaxStep = [1e10;1e10];
NumSteps = 100;
for k=1:2:length(varargin)
    a = varargin(k);
    b = varargin(k+1);
    if strcmp(a,'uinitial')
        uinitialPDE = b{1};
    elseif strcmp(a,'uo')
        uoPDE = b{1};
    elseif strcmp(a,'coeff')
        coeff = b{1};
    elseif strcmp(a,'bc')
        bc = b{1};
    elseif strcmp(a,'UpdateGeometry')
        UpdateGeometry = b{1};
    elseif strcmp(a,'Hadd')
        Hadd = b{1};
    elseif strcmp(a,'Radd')
        Radd = b{1};
    elseif strcmp(a,'Qadd')
        Qadd = b{1};
    elseif strcmp(a,'Gadd')
        Gadd = b{1};
    elseif strcmp(a,'ResidualAbsTol')
        ResidualAbsTol = b{1};
    elseif strcmp(a,'SolRelTol')
        SolRelTol = b{1};      
    elseif strcmp(a,'MaxStep')
        MaxStep = b{1};
    elseif strcmp(a,'NumSteps')
        NumSteps = b{1};        
    end
end
numNodes = size(coeff.geometry.mesh.p,2);
N = coeff.dimension;
if isempty(Hadd)
    Hadd = sparse(0,N*numNodes);
end
if isempty(Radd)
    Radd = sparse(0,1);
end
if isempty(Qadd)
    Qadd = sparse(N*numNodes,N*numNodes);
end
if isempty(Gadd)
    Gadd = zeros(N*numNodes,1);
end
% find size of H constraint matrix to determine over all dimension of
% problem
[~,~,Hmat,R] = bc.getMatrices('u',zeros(N*numNodes,1));
if isempty(find(Hmat ~= 0))
    R = [];
end
sizeH = size(R,1) + size(Radd,1);
if isempty(uoPDE)
    uoPDE = zeros(N*numNodes+sizeH,1);
end
if isempty(uinitialPDE)
    uinitialPDE = zeros(N*numNodes+sizeH,1);
    uinitialPDE(1:N*numNodes,1) = uoPDE(1:N*numNodes);
end
J = [];
    function rhs = computeRHS(u)
        uPDE = u(1:N*numNodes);
        uu = reshape(uPDE,numNodes,[]);
        p = uu(:,1:2)';
        if UpdateGeometry
            [K,M,F] = coeff.getMatrices('u',uPDE-uoPDE(1:N*numNodes),'p',p);
            [Q,G,H,R] = bc.getMatrices('u',uPDE-uoPDE(1:N*numNodes),'p',p);
        else
            [K,M,F] = coeff.getMatrices('u',uPDE-uoPDE(1:N*numNodes));
            [Q,G,H,R] = bc.getMatrices('u',uPDE-uoPDE(1:N*numNodes));
        end
        if isempty(find(H ~= 0,1))
            H = [];
            R = [];
        end
        H = [H;Hadd];
        R = [R;Radd];
        bigMatrix = [K+M+Q+Qadd H.';H sparse(sizeH,sizeH)];
        B = [F+G+Gadd;R];
        rhs = bigMatrix*(u-uo) - B;
        [jK,jM] = coeff.getJacobian('u',uPDE-uoPDE(1:N*numNodes));
        J = [jK+jM+Q+Qadd H.';H sparse(sizeH,sizeH)];
    end
uo = uoPDE;
u = uinitialPDE;
converged = false;
iter = 1;
while ~converged && iter < NumSteps
    b = computeRHS(u);
    Jxb = b;
    xb = u;
    du = J\(-b);
    maxdy = max(abs(du(numNodes+1:2*numNodes)));
    if maxdy > MaxStep(2)
        du = du*MaxStep(2)/maxdy;
    end
    maxdx = max(abs(du(1:numNodes)));
    if maxdx > MaxStep(1)
        du = du*MaxStep(1)/maxdx;
    end    
    u = u + du;
    fprintf('|du|/|u| = %e |rhs| = %e\n',norm(du)/norm(u),norm(b));
    if (norm(du)/norm(u)) < SolRelTol && norm(b) < ResidualAbsTol
        fprintf('Found geometrically consistent solution\n');
        converged = true;
    end
    iter = iter + 1;
end
end