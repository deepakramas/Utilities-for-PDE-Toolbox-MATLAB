syms sigma11 sigma22 sigma12 D1 D2;
syms eps11 eps22 gamma12 E1 E2;
syms C11 C12 e11 e31 C22 e13 e33 G12 e14 e34 Eps1 Eps2;
syms dudx dvdy dudy dvdx dphidx dphidy;

A = [C11 C12 0 e11 e31;C12 C22 0 e13 e33;0 0 G12 e14 e34;0 0 G12 e14 e34;e11 e13 e14 Eps1 0;...
    e31 e33 e34 0 Eps2];
% = [sigma11;sigma22;gamma12;D1;D2];
R = [eps11;eps22;gamma12;E1;E2];
result1 = A*R;
result2 = subs(result1,[eps11 eps22 gamma12 E1 E2],[dudx dvdy (dudy+dvdx) -dphidx -dphidy]);
E = equationsToMatrix(result2,[dudx dudy dvdx dvdy dphidx dphidy]);
EVal = subs(E,[C11 C12 C22 e11 e13 e14 e31 e33 e34 Eps1 Eps2],[C11val C12val C22val e11val e13val e14val e31val e33val e34val Eps1val Eps2val]);
function [cij11,cij12,cij21,cij22] = cijFunction()
for i=1:3
    for j=1:3
        cij11(i,k) = EVal(2*(i-1)+1, 2*(j-1)+1);
        cij12(i,k) = EVal(2*(i-1)+1, 2*(j-1)+2);
        cij21(i,k) = EVal(2*(i-1)+2, 2*(j-1)+1);
        cij22(i,k) = EVal(2*(i-1)+2, 2*(j-1)+2);
    end
end
end