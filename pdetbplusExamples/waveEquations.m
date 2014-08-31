function [equations,variables] = waveEquations(varargin)
displayEquations = true;
sigma = 0.0;
for k=1:2:length(varargin)
  a = varargin(k);
  b = varargin(k+1);
  if strcmp(a,'sigma')
      sigma = b{1};
  end
end
syms u v ddx ddy dudx dudy dvdx dvdy dudt dvdt real;
equations = [dudt - v,-ddx*dudx + sigma^2*u + sigma*2*v + dvdt];
if displayEquations
  display Problem_Equations;
  equations'
end

variables = [u v];