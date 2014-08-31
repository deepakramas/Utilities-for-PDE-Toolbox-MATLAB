function output = ts3
%% Reference
% "Quantitative studies of irradiation-induced segregation and grain boundary migration in Fe-Cr-Ni alloy"
% S. Watanabe, N. Sakaguchi, N. Hashimoto, H. Takahashi, Journal of Nuclear
% Materials 224 (1195), pp 158-168

%% domain variables
syms x y t;

%% problem variables; concentrations of three elements A,B,C and two defects v (vacancies),i (interstitials)
syms C_A(x,y,t) C_B(x,y,t) C_C(x,y,t) C_v(x,y,t) C_i(x,y,t);
Ckappa = [C_A(x,y,t);C_B(x,y,t);C_C(x,y,t)];
Cj = [C_v(x,y,t);C_i(x,y,t)];

%% total and intrinsic diffusivities

% simple hack to override built-in function "alpha()" since "syms alpha" does not
% override it
alpha = 0; 
% atomic volume, damage rate, recombination coefficient, thermodynamic
% factor, internal sink strength coefficient, boltzmann factor and
% temperature
syms Omega K_0 R alpha K_v_s K_i_s k T;
% jump distance and recombination site number
syms lambda Z_i Z_v;
Z = [Z_v,Z_i];
% migration energies
syms E_i_A(CA) E_v_A(CA) E_i_B(CB) E_v_B(CB) E_i_C(CC) E_v_C(CC);
E = [E_v_A(Ckappa(1)) E_v_B(Ckappa(2)) E_v_C(Ckappa(3));...
    E_i_A(Ckappa(1)) E_i_B(Ckappa(2)) E_i_C(Ckappa(3))];
% defect and solute pair jump frequency
syms nu_i_A nu_v_A nu_i_B nu_v_B nu_i_C nu_v_C;
syms nu_0_i_A nu_0_i_B nu_0_i_C nu_0_v_A nu_0_v_B nu_0_v_C;
Nu0 = [nu_0_v_A nu_0_v_B nu_0_v_C; nu_0_i_A nu_0_i_B nu_0_i_C];
% intrinsic diffusivities
d = sym(zeros(3,2));
for j = 1:2
    for Kappa = 1:3
        Nu = Nu0(j,Kappa)*exp(-E(j,Kappa)/(k*T));
        d(Kappa,j) = lambda^2*Z(j)*Nu/6;
    end
end
% total diffusivities
Dkappa = d*Cj;
Dj = sum(d,1)*diag(Cj);

%% flux equations 
% system 4
Jkappa = sym(zeros(2,3));
for kappa=1:3
    Jkappa(:,kappa) = 1/Omega*(-Dkappa(kappa)*gradient(Ckappa(kappa),[x y]) + ...
        d(kappa,1)*Ckappa(kappa)*gradient(Cj(1),[x,y]) - ...
        d(kappa,2)*Ckappa(kappa)*gradient(Cj(2),[x,y]));
end
Jj = sym(zeros(2,2));
for j=1:2
    if j == 1
        mySign = 1;
    else
        mySign = -1;
    end
    Jj(:,j) = 1/Omega*(mySign*(d(1,j)-d(3,j))*Cj(j)*alpha*gradient(Ckappa(1),[x,y]) + ...
        mySign*(d(2,j)-d(3,j))*Cj(j)*alpha*gradient(Ckappa(2),[x,y]) - ...
        Dj(j)*gradient(Cj(j),[x,y]));
end

%% continuity equations
% system 1
for kappa = 1:3
    eq(kappa) = -diff(Ckappa(kappa),t) - Omega*divergence(Jkappa(:,kappa),[x,y]);
end
% systems 2 and 3
syms eta C_s C_0_v C_0_i;
Ks = [K_v_s,K_i_s];
C0j = [C_0_v C_0_i];
for j = 1:2
    eq(3+j) = -diff(Cj(j),t) - divergence(Jj(:,j),[x y]) + eta*K_0 - R*prod(Cj) - Ks(j)*C_s*(Cj(j)-C0j(j));
end
output.equations = eq;
output.variables = [Ckappa;Cj];
output.internalFunctions = E(:);
end