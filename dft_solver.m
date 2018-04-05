function [box, Omega, iter, box_Picard, Omega_Picard, iter_Picard] = dft_solver( rho_guess, model, box, mesh, components, skip_Picard)
%dft_solver Summary of this function goes here
%   Detailed explanation goes here

if isnan(rho_guess)
    %initial guess for rho is just rho0*exp(-mu_ext), maxxed out at rho0
    rho_guess = min(box.rho0.*exp(-model.external.calc_mu(box, mesh, components)), box.rho0);
end
box = model.initialize_box(rho_guess,box,mesh,components);
    
% helper function for the golden-section search
function E = calc_objective_line_search(alpha)
    box_alpha = box.mix(new_box,alpha);
    E = model.calc_Picard_step(box_alpha, mesh, components);
end

min_step_multiplier = 0.0005;

err0 = model.calc_Picard_step(box, mesh, components);
err = 1;
rel_err = 1;
iter_Picard = 0;
alphas = [];
errs = [];
errs_autocorr_1 = [];
prev_step_to_new_step_errs = [];
prev_step_to_new_step_errs_autocorr_1 = [];
prev_new_vec = box.vec();
while ~skip_Picard && iter_Picard < 5000 && rel_err > 1e-5 && err > 1e-3
    % take a picard iteration step
    [err, new_box] = model.calc_Picard_step(box, mesh, components);
    
    if isinf(err0)
        rel_err = 1;
        err0  = err;
    else
        rel_err = err/err0;
    end
    
    errs = [errs, err];
    
    new_vec = new_box.vec();
    prev_step_to_new_step_errs = [prev_step_to_new_step_errs, norm(new_vec(:) - prev_new_vec(:))];
    prev_new_vec = new_vec;
    
    if iter_Picard > 11
        last10 = prev_step_to_new_step_errs(end-10:end);
        prev_step_to_new_step_errs_autocorr_1 = [prev_step_to_new_step_errs_autocorr_1, mean((last10 - mean(last10)).*(circshift(last10,1) - mean(last10)))/var(last10)];
        if prev_step_to_new_step_errs_autocorr_1(end) > 0.45
            min_step_multiplier = min_step_multiplier*1.05;
        end
    end
    
    if iter_Picard > 11
        last10 = errs(end-10:end);
        errs_autocorr_1 = [errs_autocorr_1, mean((last10 - mean(last10)).*(circshift(last10,1) - mean(last10)))/var(last10)];
        if errs_autocorr_1(end) < 0
            min_step_multiplier = min_step_multiplier*0.95;
        end
    end
    
    % enforce a minimum min_step_multiplier
    if min_step_multiplier < 0.0001
        min_step_multiplier = 0.0001;
    end
    
    % line search in the direction of new_box-old_box
    old_infnorm_n3 = norm(box.n3(:),Inf);
    new_infnorm_n3 = norm(new_box.n3(:),Inf);
    
    % the largest step we can take would make the packing fraction (n3)
    % equal 0.9 (since n3 >= about 1 would be unphysical)
    alpha_max = min(0.9,(0.9-old_infnorm_n3)/(new_infnorm_n3-old_infnorm_n3));
    if alpha_max < 0
        alpha_max = 0.9;
    end
    
    % bracket the line search, to some extent
    err_alpha = calc_objective_line_search(alpha_max);
    while err_alpha > 100*err || isnan(err_alpha)
        alpha_max = alpha_max/2;
        err_alpha = calc_objective_line_search(alpha_max);
    end
        
    % line search
    a = 0;
    b = alpha_max;
    min_step = b*min_step_multiplier;
%     a = golden_section_search( @calc_objective_gss, a, b, tol_gss );
    a = three_point_quadratic_search( @calc_objective_line_search, a, b, err, err_alpha, min_step );
    alphas = [alphas a];
    box = box.mix(new_box,a);

    iter_Picard=iter_Picard+1;
    % error display
    if mod(iter_Picard,100) == 0 || iter_Picard == 1
        fprintf('iter = %4d, err = %9.4g, rel_err = %9.4g, min_step_multiplier = %9.4g\n', iter_Picard, err, rel_err, min_step_multiplier)
        
        figure(2)
        subplot(2,3,2)
        plot(prev_step_to_new_step_errs_autocorr_1)
%         plot(mesh.x,box.rho(:,1,1,1),'ro-',mesh.x,box.rho(:,1,1,2),'bo-',mesh.x,new_box.rho(:,1,1,1),'r--',mesh.x,new_box.rho(:,1,1,2),'b--')
        subplot(2,3,4)
        semilogy(1:iter_Picard,errs,'k',1:iter_Picard,prev_step_to_new_step_errs,'g')
%         plot(mesh.x,box.D_field(:,1,1,1),'ro-',mesh.x,box.D_field(:,1,1,2),'bo-')
        subplot(2,3,3)
        plot(errs_autocorr_1)
%         rho_seg = exp(box.mu_bath(:,:,:,components.seg_to_poly) + box.P_f + box.P_b + box.D_field(:,:,:,components.seg_to_type) - box.V_external(:,:,:,components.seg_to_type));
%         plot(mesh.x(:,1,1),squeeze(rho_seg(:,1,1,:)))
        subplot(2,3,1)
        plot(alphas)
        subplot(2,3,5)
        if mesh.D == 1
            plot(mesh.x,squeeze(box.rho))
        elseif mesh.D == 2
            imagesc(box.rho(:,:,1,1))
        end
        drawnow
    end
end
if mod(iter_Picard,100) ~= 0 && iter_Picard > 1
        fprintf('iter = %4d, err = %9.4g, rel_err = %9.4g, min_step_multiplier = %9.4g\n', iter_Picard, err, rel_err, min_step_multiplier)
end

Omega_Picard = model.calc_Omega(box, mesh, components)
box_Picard = box;

% helper function for the KNL optimization
function [F] = calc_F(x)
    box = box.unvec(x,mesh,components);
    new_box = model.calc_Newton_step(box, mesh, components);
    F = new_box.vec() - x;
end

guess = box.vec();
[sol, iter, ierr] = knl(guess,@calc_F,knl_optset('debug',1,'maxit',50,'rtol',1e-10,'orthog','cgs'));
if ierr > 0
    if ierr==1 % hit maximum iterations
        knl_err = iter(end,1)
    else % ierr==2 line search fail
        knl_err = iter(end-1,1)
    end
    if knl_err > 1e-4
        error('knl failure: no good solution reached')
    end
    warning('knl failure')
end
box = box.unvec(sol,mesh,components);

% [Omega, Omega_ideal_gas, Omega_external, Omega_mu, Omega_excess, Omega_polymer] = model.calc_Omega(box, mesh, components)
Omega = model.calc_Omega(box, mesh, components)

end

