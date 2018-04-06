clear

% This is an example of a pure hard sphere fluid in a 2D square pore
% surrounded by hard walls, with parameters chosen to match Figure 2 of
% Roth, J. Phys.: Condens. Matter 22 (2010) 063102.

% 2D symmetry, periodic BCs, cell size of 1/20 sigma, 400 by 400 cells
% (i.e. the box is 20 by 20 sigma)
mesh = dft_mesh_parameters(2, 'p', 0.05, 400,400)

% one component, HS diameter of 1 sigma, packing fraction of 0.3665
% converted to number density
components = dft_component_parameters(1, 1, 0.3665*6/pi)

% Box initialization
% the walls exist along x=0 and y=0, which due to periodic BCs, also means
% they are at x=20 and y=20
wall_location = logical(zeros(400,400)); %#ok<LOGL>
wall_location(1,:) = true;
wall_location(:,1) = true;
box = dft_box(mesh, components, wall_location)

% external potential for the hard wall
vext = zeros(400,400);
vext(1:10,:) = 200;
vext(:,1:10) = 200;
vext(392:400,:) = 200;
vext(:,392:400) = 200;
external = dft_functional_external(vext, mesh, components)

% ideal gas and hard sphere functionals
id = dft_functional_ideal_gas
hs = dft_functional_hard_spheres_WhiteBear

% basic model to put everything together
model = dft_model_basic(id, external, hs,mesh,components)

% default initial guess
rho_guess = NaN;

% call the solver, skipping picard iteration
[box, Omega, iter] = dft_solver( rho_guess, model, box, mesh, components, true )

% plot the results as in Figure 2 of Roth 2010
figure(1)
diag_rho = diag(box.rho);
subplot(1,2,1), plot(sqrt(2)*0.05*(0:1:40)-sqrt(2)/2+0.5,diag_rho(1:41))
diag2_rho = diag(box.rho(15:200,:,1,1));
subplot(1,2,2), plot(sqrt(2)*0.05*(0:1:40)-sqrt(2)/2+0.5,diag2_rho(1:41))
