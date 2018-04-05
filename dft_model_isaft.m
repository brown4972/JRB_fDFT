classdef dft_model_isaft < dft_model
    %dft_model_isaft - class to store functionals and calculate Omega and mu
    %   This is version of the DFT equations as written for iSAFT
    %   see Jain et al. 2007
    
    properties (SetAccess = protected)
        polymer = dft_functional; % polymer functional
    end
    
    methods
        function obj = dft_model_isaft(ideal_gas, external, excess, polymer, mesh, components)
            obj = obj@dft_model(ideal_gas, external, excess, mesh, components);
            
            if ~isa(polymer,'dft_functional_isaft')
                error('ERROR: polymer must be a dft_functional_isaft')
            end
            obj.polymer = polymer;
        end
        
        function box = initialize_box(obj, rho_guess, box, mesh, components)
            zero_Ncomp = zeros(mesh.NX,mesh.NY,mesh.NZ,components.Ncomp);
            zero_Npoly = zeros(mesh.NX,mesh.NY,mesh.NZ,components.Npoly);
            
            box_bath = box.set_rho(box.rho0,mesh,components);
            box_bath = box_bath.set_static_fields(zero_Npoly, zero_Ncomp);
            
            mu_bath_excess = arrayfun(@(x)(x.calc_mu(box_bath,mesh,components)),obj.excess,'UniformOutput',false);
            mu_bath_excess = cat(5,mu_bath_excess{:}); % creates an NX by NY by NZ by N by length(excess) array that stores the mu's of the various excess functionals
            mu_bath_excess = sum(mu_bath_excess,5);
                                    
            mu_bath = obj.polymer.calc_mu_bath(box_bath, mu_bath_excess, mesh,components);
            V_external = obj.external.calc_mu(box_bath, mesh, components);
            
            
            box = initialize_box@dft_model(obj, rho_guess, box, mesh, components);
            box = box.set_static_fields(mu_bath, V_external);
        end
        
        function [Omega, Omega_ideal_gas, Omega_external, Omega_mu, Omega_excess, Omega_polymer] = calc_Omega(obj, box, mesh, components)
            [Omega, Omega_ideal_gas, Omega_external, Omega_mu, Omega_excess] = calc_Omega@dft_model(obj, box, mesh, components);

            mu_excess = arrayfun(@(x)(x.calc_mu(box,mesh,components)),obj.excess,'UniformOutput',false);
            mu_excess = cat(5,mu_excess{:}); % creates an NX by NY by NZ by Ncomp by length(excess) array that stores the mu's of the various excess functionals
            mu_excess = sum(mu_excess,5);
            Omega_polymer = obj.polymer.calc_free_energy(mu_excess, box, mesh, components);
            
            Omega = Omega + Omega_polymer;
        end
        
        % Omega_mu is the contribution to the grand free energy from the
        % fixed bath chemical potential
        function [ Omega_mu ] = calc_Omega_mu(obj, box, mesh, components)
            Omega_mu = -sum(box.weighted_rho_seg.*box.mu_bath(:,:,:,components.seg_to_poly),4);
            Omega_mu = Omega_mu.*mesh.cell_weight*mesh.DX^(mesh.D);
            Omega_mu = sum(Omega_mu(:));
        end
                    
        function [err, new_box] = calc_Picard_step(obj, box, mesh, components)    
            % for Picard iteration, we regenerate the fields at each step
            % from the current density profiles...
            mu_excess = arrayfun(@(x)(x.calc_mu(box,mesh,components)),obj.excess,'UniformOutput',false);
            mu_excess = cat(5,mu_excess{:}); % creates an NX by NY by NZ by Ncomp by length(excess) array that stores the mu's of the various excess functionals
            mu_excess = sum(mu_excess,5);

            [D_field, P_f, P_b] = obj.polymer.calc_fields(mu_excess,box,mesh,components);

            % and use the resulting fields to generate a new set of density
            % profiles
            new_rho_seg = exp(box.mu_bath(:,:,:,components.seg_to_poly) + P_f + P_b ...
                                + D_field(:,:,:,components.seg_to_type) - box.V_external(:,:,:,components.seg_to_type));
            
            % This effectively makes the decision variable log(rho) instead
            % of rho in the Picard iteration. This is useful since log(rho)
            % is unbounded, but rho must be non-negative
            err = log(box.rho_seg) - log(new_rho_seg);
            err = norm(err(:));

            if nargout > 1
                new_box = box.set_rho(new_rho_seg,mesh,components);
            end
        end
        
    end
    
end

