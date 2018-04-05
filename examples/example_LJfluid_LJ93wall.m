clear

% This is an example of an attractive Lennard-Jones (LJ) fluid near a
% attractive LJ (soft) wall. This example also uses reflective instead of
% periodic BCs

% When using reflective BCs instead of periodic, the cells on the edges of
% the box are NOT replicated/reflected like all the cells inside
% Equivalently, the cells on the edges are "half cells" (on the corners in
% 2D "quarter cells", etc) whena re then reflected to make full cells.

% 1D symmetry, reflective BCs, cell size of 1/40 sigma, box size of 1001
% cells (box size is (1/2 + 999 + 1/2) * 1/40 = 25 sigma) 
mesh = dft_mesh_parameters(1, 'r', 0.025, 1001)

% one component, HS diameter 1 sigma, bath density of 0.846/sigma^3
components = dft_component_parameters(1, 1, 0.846)

% Set up the box
% There is no wall weight here because there is no hard wall
box = dft_box(mesh, components)

% Set the external potential to the the 9-3 Lennard-Jones potential,
% maqxxed out at 200kT
% the parameters for the LJ potential are:
% epsilon=2.1kT, sigma_wall=1, cut off at 3 sigma
vext = min(u_LJ93_CS(mesh.x,[2.1,1,3]),200);
external = dft_functional_external(vext, mesh, components)

% ideal gas and hard sphere functionals
id = dft_functional_ideal_gas
hs = dft_functional_hard_spheres_FMT

% This is the attractive LJ functional; the first argument is a funtion
% handle to the pairwise potential, and the second arguement to the
% functioal class is a cell array of matrices that set the parameters of
% the choosen pairwise potential
epsilon = 1;
sigma = 1;
r_cut = 2.5;
lj = dft_functional_pairwise_generic(@u_LJ_WCA, {epsilon,sigma,r_cut}, mesh, components)

% Note that the third argument to the model is an array of the various
% excess free energy functionals (since there can be any number of these)
model = dft_model_basic(id, external, [hs,lj],mesh,components)

% default initial guess
rho_guess = NaN;

% call the solver w/o skipping Picard
[box, Omega, iter] = dft_solver( rho_guess, model, box, mesh, components, false )

% plot the answer
figure(1)
plot(mesh.x,box.rho,'.-')
