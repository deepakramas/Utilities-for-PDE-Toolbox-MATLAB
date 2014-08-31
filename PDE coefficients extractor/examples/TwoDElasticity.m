% Define equations for geometrically nonlinear linearly elastic plane stress material
function output = TwoDElasticity
%% Modulus of Elasticity, Poisson's ration, Density, constant forcing terms
syms El nu rho f1 f2;
output.parameters = [El,nu,rho,f1,f2];
syms u(x,y,t) v(x,y,t) real;
u = u(x,y,t);
v = v(x,y,t);
output.variables = [u,v];
%% material constants
G = El/(1-nu^2);
C11 = G;
C12 = G*nu;
C22 = G;
C66 = G*(1-nu)/2;
%% deformation gradient
F = sym(eye(2)) + jacobian([u v],[x y]);%[dudx dudy;dvdx dvdy];
%% Green-Lagrange strain tensor
epsilon =  1/2*(F.'*F - sym(eye(2)));
epsilon11 = epsilon(1,1);
epsilon22 = epsilon(2,2);
epsilon12 = epsilon(1,2);
%% 2nd Piola–Kirchhoff stress tensor
Sx = C11*epsilon11 + C12*epsilon22;
Sy = C12*epsilon11 + C22*epsilon22;
Txy = C66*2*epsilon12; % watch out for factor of 2 for engineering strain
PK2 = [Sx Txy;Txy Sy];
PK1 = PK2*F.'; 
output.divergenceTerm = -PK1;
output.equations = [divergence(-PK1(:,1),[x y]), divergence(-PK1(:,2),[x y])] + ...
    rho*[diff(u,t,t), diff(v,t,t)] - [f1,f2];
output.symmetricC = false;
end