clear

% This is an example of a pure hard sphere fluid next to a hard wall, with
% parameters chosen to match Figure 1(a) of Roth, J. Phys.: Condens. Matter
% 22 (2010) 063102.

% The mesh parameters class stores static properties of the computational
% mesh, including it's size and boundary conditions
% In this case, the system has 1D symmetry, the bounary conditions are
% periodic, each cell is 1/40th of sigma (the nondimensionalized unit of
% distance), and the box is 10 sigma wide
dimensionality = 1
boundary_conditions = 'p' %periodic
cell_size = 0.025
number_of_cells = 400
mesh = dft_mesh_parameters(dimensionality, boundary_conditions, cell_size, number_of_cells)


% The component parameters class stores static properties of the various
% hard sphere components in the system.
% In this cae we have a single component of diameter 1 sigma that's in
% equilibrium with a "bath" of number density 0.813/sigma^3
Num_components = 1
hard_sphere_diameters = [1]
rho_bulk = 0.813
components = dft_component_parameters(Num_components, hard_sphere_diameters, rho_bulk)


% The box class stores the non-static solution to the fDFT equations
% The location of the wall is passed in here
% Note that the wall is at cell 1, which due to periodic BCs is both x=0
% and x=10
wall_location = logical(zeros(number_of_cells,1)); %#ok<LOGL>
wall_location(1) = true;
box = dft_box(mesh, components, wall_location)


% The external potential functional
% in this case the wall at x=0 and x=10
vext = zeros(number_of_cells,1);
vext(1:20) = 200;
vext(382:number_of_cells) = 200;
external = dft_functional_external(vext, mesh, components)

% The ideal gas functional has no input parameters
id = dft_functional_ideal_gas

% There are several options for hard sphere functionals, this one is 
% Roth's original fundamental measure theory
hs = dft_functional_hard_spheres_FMT

% The model class puts the various functionals together
model = dft_model_basic(id, external, hs, mesh, components)

% This is the initial guess for the density profiles.
% Setting the initial guess to NaN means that it will be set to
% rho_bulk*exp(-v_ext) 
rho_guess = NaN;

% This is the solver call, the only parameter to set if if we want to skip
% the (slow) Picard solver, if we think our inital guess is close enough
% not to need it
skip_Picard = true;
[box, Omega, iter] = dft_solver( rho_guess, model, box, mesh, components, skip_Picard )

% plot the solution
figure(1)
plot(mesh.x,box.rho,'.-')
