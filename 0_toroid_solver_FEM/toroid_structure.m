% toroid parameter

function y=toroid_structure

toroid_geom.lambda=1.55e-6;
toroid_geom.n_guide=1.4457; % 1.4457 1.4447453
toroid_geom.n_sub=1.0; % air@1.0  water@1.33
toroid_geom.n_metal=0.14447+11.366i; % silver@0.63 0.0245+4.7635i; silver@1.55 0.14447+11.366i
toroid_geom.core_radius=2e-6;
toroid_geom.toroid_R=30e-6; %3e-5
toroid_geom.metal_thickness=0.1e-6; %10~100 nm
toroid_geom.metal_radius=toroid_geom.core_radius+toroid_geom.metal_thickness;
toroid_geom.window=3.0*toroid_geom.core_radius;
toroid_geom.num_modes=10;%round(1.2*V_num^2/2);
toroid_geom.k0_a=(2*pi*toroid_geom.core_radius/toroid_geom.lambda);
toroid_geom.M_guess=177;

y=toroid_geom;

end