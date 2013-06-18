% toroid FEM mode slover based on Oxborrrow's code
% 2013.06

clear all;
close all;

c_m_per_sec=299792458;
u0=4*pi*1e-7;
e0=8.8542e-12;

% toroid parameter

toroid_geom=toroid_structure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fem_tor=toroid_mode_solver_weak(toroid_geom); % toroid mode solving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sol_frequency,sol_ind]=sort(fem_tor.sol.lambda,'ascend');
sol_wavelength=c_m_per_sec./sol_frequency.*1e9;

