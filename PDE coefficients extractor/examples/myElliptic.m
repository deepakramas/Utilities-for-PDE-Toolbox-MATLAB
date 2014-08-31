function output = myElliptic
syms u(x,y)
u = u(x,y);
output.variables = u;
output.equations = divergence(u*(sqrt(x)/sqrt(u)+sqrt(u)*4*x/y)*gradient(u,[x y]),[x y]) + u - 3;
end