clear

% This is an example of a 10-site symmetric diblock copolymer forming a
% lamellar morphology. This shows an example of setting a specific initial
% guess.

% 1D symmetry, periodic BCs, cell spacing of 1/10 sigma, 50 cell box size
% (i.e. 5 sigma)
mesh = dft_mesh_parameters(1, 'p', 0.1, 50)

% 2 monomer types, HS diameter of 1 and 1, 1 polymer type, polymer bath
% density of 0.85/sigma^3, polymer sequence is 5 of type 1 followed by 5 of
% type 2
components = dft_polymer_parameters(2, [1 1], 1, 0.85, {[ 1 1 1 1 1 2 2 2 2 2 ]})
box = dft_box_isaft(mesh, components)

% There is no external potential, so this is set to the dummy functional
external = dft_functional

% ideal gas and hard sphere functionals
id = dft_functional_ideal_gas
hs = dft_functional_hard_spheres_asFMT

% Attractive LJ interactions. The like-like interactions are stronger than
% the unlike interactions, this drives the phase separation. Note that the
% parameters are passed in as symmetric matrices.
epsilon = [1 0.8; 0.8, 1];
sigma = ones(2);
r_cut = 2.5*ones(2);
dHS = ones(2);
lj = dft_functional_pairwise_generic(@u_LJ_WCA, {epsilon,sigma,r_cut}, mesh, components)

% iSAFT polymer funtional
polymer = dft_functional_isaft(components)

% iSAFT polymer model
model = dft_model_isaft(id, external, [hs, lj], polymer, mesh, components)

% Initial guess is an array with dimensions X size by Y size by Z size by
% number of monomer types. The box is initialized with the left half rich
% in type 2, the right side rich in type 1, and the total density equal to
% the bath density
rho_guess = 0.8*ones(50,1,1,2);
rho_guess(1:25,1,1,1) = 0.05;
rho_guess(26:50,1,1,2) = 0.05;

% call the solver with Picard
[box, Omega, iter] = dft_solver( rho_guess, model, box, mesh, components, false )

% plot the results
figure(1)
plot(mesh.x,box.rho(:,1,1,1),'ro-',mesh.x,box.rho(:,1,1,2),'bo-')
