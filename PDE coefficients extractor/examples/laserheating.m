function output = laserheating
syms E(z,r,t)
syms phi
syms omega omega_e omega_0
syms rho
syms k K beta_ n_2 c 
E = E(z,r,t);
%rho = rho(t);
%omega = omega(rho);
%omega_epsilon = omega_epsilon(rho);
i = complex(0,1);
output.equations = i/(2*k)*diff(E,r,r)+i*omega/c*n_2*abs(E)^2*E-i*k*omega_e^2/omega_0^2*E-beta_^K/2*abs(E)^(2*K-2)*E-diff(E,z);
output.coordinate.cylindrical.radial = r;
output.coordinate.cylindrical.azimuth = phi;
output.coordinate.cylindrical.height = z;
output.variables = E;
output.parameters = [k K beta_ n_2 c rho omega omega_e omega_0];
end