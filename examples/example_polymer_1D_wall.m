clear

% 20 site polymer between two walls

% 1D symmetry, periodic BCs, cell size of 1/40 sigma, box size of 220 cells (5.5 sigma)
mesh = dft_mesh_parameters(1, 'p', 0.025, 220)

% For polymers the component creation is a bit different, because we now
% can have different types of monomers, and different types of polymers
% (that may have sequences containing multiple types of monomers)
% Note that monomers don't exist on their own, so if we are using the
% polymer functionals, lone hard sphere species need to be "polymers" of
% length one
% In this case, we have a single monomer type and a single polymer type
Num_monomers = 1;
monomer_dHS = 1;
Num_polymers = 1;
rho_polymer_bath = 0.813;
polymer_sequences = {ones(1,20)};
components = dft_polymer_parameters(Num_monomers, monomer_dHS, Num_polymers, rho_polymer_bath, polymer_sequences)

% wall at x=0 aka x=5.5
wall_location = logical(zeros(220,1)); %#ok<LOGL>
wall_location(1) = true;
box = dft_box_isaft(mesh, components, wall_location)

% external field for wall(s)
vext = zeros(220,1,1,1);
vext(1:20,:,:,:) = 200;
vext(202:220,:,:,:) = 200;
external = dft_functional_external(vext, mesh, components)

% ideal gas and hard sphere functionals are unchanged
id = dft_functional_ideal_gas
hs = dft_functional_hard_spheres_asFMT

% polymer functional needs to know about the sequences etc
polymer = dft_functional_isaft(components)

% use the isaft model instead of the basic model, since isaft rearranges
% the DFT equations a bit
model = dft_model_isaft(id, external, hs, polymer, mesh, components)

% default initial guess
rho_guess = NaN;

% call the solver without skipping the Picard iteration
[box, Omega, iter] = dft_solver( rho_guess, model, box, mesh, components, false )

% plot the results
figure(1)
plot(mesh.x,box.rho)
