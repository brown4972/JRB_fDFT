clear

% This is an example of a binary mixture of two different sized hard sphere
% components (size 1 and 3 sigma) between two hard walls, with parameters
% chosen to match Figure 3 of Roth, J. Phys.: Condens. Matter 22 (2010)
% 063102.

% 1D symmetry, periodic BCs, cell size of 1/40 sigma, 2000 cells (50 sigma)
mesh = dft_mesh_parameters(1, 'p', 0.025, 2000)

% 2 components, HS diameters of 1 and 3, bath densities as specified to
% give packing fractions of 0.01 and 0.3576, respectively.
components = dft_component_parameters(2, [1,3], [0.01*6/pi, 0.3576*6/pi/27])

% wall at x=0 and x=50
wall_location = logical(zeros(2000,1)); %#ok<LOGL>
wall_location(1) = true;
box = dft_box(mesh, components, wall_location)

% external potentials for walls at x=0 and x=50, depend on the size of the
% components
vext = zeros(2000,1,1,2);
vext(1:20,1,1,1) = 200;
vext(1982:2000,1,1,1) = 200;
vext(1:60,1,1,2) = 200;
vext(1942:2000,1,1,2) = 200;
external = dft_functional_external(vext, mesh, components)

% ideal gas and hard sphere functionals
id = dft_functional_ideal_gas
hs = dft_functional_hard_spheres_asFMT

% model to put functionals together
model = dft_model_basic(id, external, hs, mesh, components)

% default initial guess
rho_guess = NaN;

% call solver without skipping picard iteration
[box, Omega, iter] = dft_solver( rho_guess, model, box, mesh, components, false )

% plot results
figure(1)
subplot(1,2,1), plot(mesh.x,box.rho(:,1,1,1))
xlim([0 8])
subplot(1,2,2), plot(mesh.x,box.rho(:,1,1,2))
xlim([0 8])
