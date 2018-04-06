clear

% This is and example of an asymmetric diblock copolymer that phase
% separates into cylinders. Unlike the lamellae example, this example uses
% a square repulsion to drive phase separation.

% 2D, periodic, cells are 1/4 sigma, 40 by 70 cells (10 by 17.5 sigma)
mesh = dft_mesh_parameters(2, 'p', 0.25, 40,70)

% 2 monomer types, both size 1, 1 polymer type, bath density 0.75/sigma^3,
% polymer is 5 type 1 then 15 type 2
components = dft_polymer_parameters(2, [1 1], 1, 0.75, {[ones(1,5), 2*ones(1,15)]})

% initialize box
box = dft_box_isaft(mesh, components)

% no external field
external = dft_functional

% ideal gas and HS functionals
id = dft_functional_ideal_gas
hs = dft_functional_hard_spheres_asFMT

% Pairwise interactions are a square resulsions (i.e. a square well with a
% negative well depth).
epsilon = [0 -0.1; -0.1, 0];
r_cut = 2*ones(2);
sw = dft_functional_pairwise_SW(epsilon,r_cut, mesh, components)

% iSAFT polymer functional
polymer = dft_functional_isaft(components)

% iSAFT model
model = dft_model_isaft(id, external, [hs, sw], polymer, mesh, components)

% initial guess: type 1 is dense in the corners and in the middle of the
% box
rho_guess = zeros(40,70,1,2);
rho_A = 0.01*ones(40,70,1,1);
rho_A(mesh.r<2.6) = 0.74;
rho_A(sqrt((mesh.x-5).^2+(mesh.y-8.75).^2)<2.6) = 0.74;
rho_guess(:,:,1,1) = rho_A;
rho_guess(:,:,1,2) = 0.75 - rho_A;

% call the solver skipping Picard iteration
[box, Omega, iter] = dft_solver( rho_guess, model, box, mesh, components, true )

% plot the results
figure(1)
surf(box.rho(:,:,1,1))
