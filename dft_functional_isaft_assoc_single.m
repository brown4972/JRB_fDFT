classdef dft_functional_isaft_assoc_single < dft_functional
    %dft_functional_isaft - iSAFT association functional
    %   this is a simplified version of iSAFT for associating fluids (i.e.
    %   breakable bonds) where each component (i.e. segment type in polymer
    %   systems) has at most only a single association site. A more general
    %   version of this functional is described in Bymaster and Chapman,
    %   JPCB (2010).
    
    
    properties
        assoc_with_type = []; % length Ncomp array mapping associations (types that have no association site are set to 0)
        prefactor = []; % prefactor that depends on association bond energies and a geometrical constant based on r_c and theta_c; length Ncomp array
    end
    
    methods
        function obj = dft_functional_isaft_assoc_single(assoc_with_type, eps_assoc, r_c, theta_c, components)
            if ~isa(components,'dft_polymer_parameters')
                error('ERROR: components must be a dft_polymer_parameters')
            end
            
            if length(assoc_with_type) ~= components.Ncomp || length(eps_assoc) ~= components.Ncomp || length(r_c) ~= components.Ncomp || length(theta_c) ~= components.Ncomp
                error('ERROR: assoc_with_type, eps_assoc, r_c, and theta_c must have length Ncomp')
            end
            
            for i=1:components.Ncomp
                j = assoc_with_type(i);
                if j ~= 0
                    if eps_assoc(i) ~= eps_assoc(j)
                        error('ERROR: eps_assoc must be symmetric wrt associations')
                    elseif r_c(i) ~= r_c(j)
                        error('ERROR: r_c must be symmetric wrt associations')
                    elseif theta_c(i) ~= theta_c(j)
                        error('ERROR: theta_c must be symmetric wrt associations')
                    end
                end
            end
            
            obj.assoc_with_type = assoc_with_type;
            
            dHS = zeros(components.Ncomp,components.Ncomp);
            for i=1:components.Ncomp
                for j=1:components.Ncomp
                    dHS(i,j) = (components.dHS(i)+components.dHS(j))/2;
                end
            end
            
            obj.prefactor = zeros(1,components.Ncomp);
            for i=1:components.Ncomp
                j = assoc_with_type(i);
                if j ~= 0
                    energy_term = (exp(eps_assoc(i))-1);
                    rc_minus_dHS = r_c(i) - dHS(i, j);
                    kappa = ( 1 - cos(theta_c(i)) )^2/4;
                    obj.prefactor(i) = energy_term*rc_minus_dHS*kappa;
                end
            end
        end
        
        function F = calc_free_energy(obj, box, mesh, components)
            sqrt_y_cav = obj.calc_y_cav_terms(box, mesh, components);
            unbonded_fraction = obj.calc_unbonded_fraction(sqrt_y_cav, box, mesh, components);
            
            % this is equation (16) from Bymaster & Chapman, simplified to
            % single association sites
            F = sum(box.weighted_rho.*(log(unbonded_fraction) - unbonded_fraction/2 + 0.5),4); % note that the contribution of unassociating components is 0
            F = F.*mesh.cell_weight*mesh.DX^(mesh.D);
            F = sum(F(:));
        end
        
        function mu = calc_mu(obj,box,mesh,components)
            [sqrt_y_cav, inv_y_cav, dy_dxi2, dy_dxi3] = obj.calc_y_cav_terms(box, mesh, components);

            unbonded_fraction = obj.calc_unbonded_fraction(sqrt_y_cav, box, mesh, components);
            
            % FIXME: vectorize?
            % this is equation (25) from Bymaster & Chapman, simplified to
            % single association sites (note my i is their j and visa versa)
            mu = log(unbonded_fraction); % note that this sets mu of unassociating components to 0
            for i=1:components.Ncomp
                for j=1:components.Ncomp
                    k = obj.assoc_with_type(j);
                    if k ~= 0
                        y_cav_xi2_term =  box.weighted_rho(:,:,:,j).*(1-unbonded_fraction(:,:,:,j)).*inv_y_cav(:,:,:,j,k).*dy_dxi2(:,:,:,j,k);
                        y_cav_xi3_term =  box.weighted_rho(:,:,:,j).*(1-unbonded_fraction(:,:,:,j)).*inv_y_cav(:,:,:,j,k).*dy_dxi3(:,:,:,j,k);
                        
                        y_cav_xi2_term = box.WXI2(:,:,:,i).*mesh.forward_fft(y_cav_xi2_term);
                        y_cav_xi3_term = box.WXI3(:,:,:,i).*mesh.forward_fft(y_cav_xi3_term);
                        
                        mu(:,:,:,i) = mu(:,:,:,i) - 0.5*mesh.inverse_fft(y_cav_xi2_term + y_cav_xi3_term);
                    end
                end
            end            
        end
        
    end
    
    methods (Access = protected)
        function [sqrt_y_cav, inv_y_cav, dy_dxi2, dy_dxi3] = calc_y_cav_terms(obj, box, mesh, components)
            Delta = 1-box.xi3;
            
            const = zeros(1,1,1,components.Ncomp,components.Ncomp);
            for i=1:components.Ncomp
                for j=1:components.Ncomp
                    d_i = components.dHS(i);
                    d_j = components.dHS(j);
                    const(1,1,1,i,j) = d_i*d_j/(d_i+d_j);
                end
            end
            
            C_xi2_over_Delta = const.*box.xi2./Delta;
            y_cav = (1 + 3*C_xi2_over_Delta + 2*C_xi2_over_Delta.^2)./Delta;
            
            sqrt_y_cav = sqrt(y_cav);
            
            if nargout > 1
                inv_y_cav = 1./(y_cav);
                
                dy_dxi2 = const.*(3+4*C_xi2_over_Delta)./Delta.^2;
                dy_dxi3 = (1+6*(C_xi2_over_Delta + C_xi2_over_Delta.^2))./Delta.^2;
            end
        end
        
        % FIXME: ?? move this function into dft_box_isaft and the input
        % parameters of this class into dft_polymer_parameters, so that the
        % unbonded_fraction can be saved along with rho, and recalculated
        % from the last iteration as a starting point
        function unbonded_fraction = calc_unbonded_fraction(obj, sqrt_y_cav, box, mesh, components)
            % initialize unbonded_fraction as all ones
            unbonded_fraction = ones(mesh.NX, mesh.NY, mesh.NZ, components.Ncomp);

            % helper function for the nonlinear solver
            function resid = calc_resid_unbonded_fraction(X)
                X = reshape(exp(X),mesh.NX, mesh.NY, mesh.NZ, components.Ncomp);
                resid = zeros(mesh.NX, mesh.NY, mesh.NZ, components.Ncomp);
                for i=1:components.Ncomp
                    j = obj.assoc_with_type(i);
                    if j ~= 0
                        integrand = box.weighted_rho(:,:,:,j).*X(:,:,:,j).*sqrt_y_cav(:,:,:,i,j);
                        integrand = mesh.forward_fft(integrand).*box.Wassoc(:,:,:,i,j);
                        integral = obj.prefactor(i).*sqrt_y_cav(:,:,:,i,j).*mesh.inverse_fft(integrand);
                        resid(:,:,:,i) = X(:,:,:,i).*(1+integral)-1;
                    else
                        resid(:,:,:,i) = 0;
                    end
                end
                resid = resid(:);
            end
            
            [soln, iter, ierr] = knl(log(unbonded_fraction(:)),@calc_resid_unbonded_fraction,knl_optset());
            
            if ierr > 0
               disp(ierr)
               warning('knl failed to solve for unbonded_fraction')
            end
            
%             soln = unbonded_fraction(:);
%             resid = calc_resid_unbonded_fraction(soln);
%             err0 = norm(resid);
%             err = 1;
%             soln = soln - resid;
%             iter = 0;
%             while err > 1e-6 && iter < 100
%                 resid = calc_resid_unbonded_fraction(soln);
%                 err = norm(resid)/err0;
%                 soln = soln - resid;
%                 iter = iter+1;
%             end
            
            unbonded_fraction = reshape(exp(soln),mesh.NX, mesh.NY, mesh.NZ, components.Ncomp);
            
        end
    end
    
end

