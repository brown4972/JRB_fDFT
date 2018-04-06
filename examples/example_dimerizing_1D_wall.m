clear

% This is an example of using the association functional (i.e. "breakable
% bonds") from Bymaster and Chapman JPCB (2010) to model a system where
% each bead has a single association site, meaning that the beads can come
% together to form dimers. The parameters chosen here match the solid blue
% curve in Figure 2(a) of that paper.

% 1D symmetry, reflective BCs, cell size of 1/40 sigma, box size of 201
% cells, equal to (1/2 + 199 + 1/2)*1/40 = 5 sigma
mesh = dft_mesh_parameters(1, 'r', 0.025, 201)

% 1 monomer type, HS size of 1, 1 polymer type, density of 0.1999/sigma^3,
% polymers consist of a single monomer
components = dft_polymer_parameters(1, 1, 1, 0.1999, {1})

% box initialization with a wall at x=0
wall_location = logical(zeros(201,1)); %#ok<LOGL>
wall_location(1) = true;
box = dft_box_isaft(mesh, components, wall_location)

% external field functional
vext = zeros(201,1,1,1);
vext(1:20,:,:,:) = 200;
external = dft_functional_external(vext, mesh, components)

% ideal gas and hard sphere functionals
id = dft_functional_ideal_gas
hs = dft_functional_hard_spheres_asFMT

% polymer functional
% In this case, since all "polymers" are single monomers, this isn't going
% to do much, but it is required because we are using the assocation
% functional also.
polymer = dft_functional_isaft(components)

% association functional
% The association functional used is a simplified one that only allows a
% single association site per type, so we just need to list what associates
% with what, the energy, and geometric constants.
assoc_with_type = [1]; %#ok<NBRAK>
eps_assoc = 14;
r_cut = 1.05; % this is the size of the association zone
theta_cut = 27*pi/180; % this is the solid angle of the association zone
assoc = dft_functional_isaft_assoc_single(assoc_with_type, eps_assoc, r_cut, theta_cut, components)

% The iSAFT model is used to put everything together, note that the hard
% sphere and association functionals are in the excess free energy list.
model = dft_model_isaft(id, external, [hs, assoc], polymer, mesh, components)

% default initial guess
rho_guess = NaN;

% call the solver, skipping picard iteration (the initial guess is actually
% pretty close to the solution)
[box, Omega, iter] = dft_solver( rho_guess, model, box, mesh, components, true )

% plot the results
figure(1)
plot(mesh.x,box.rho)
