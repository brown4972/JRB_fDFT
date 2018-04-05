classdef dft_model_basic < dft_model
    %dft_model - class to store functionals and calculate Omega and Picard
    %update steps
    %   This is the basic DFT model without polymers
    
    properties
    end
    
    methods
        function obj = dft_model_basic(ideal_gas, external, excess, mesh, components)
            obj = obj@dft_model(ideal_gas, external, excess, mesh, components);
        end
        
        function box = initialize_box(obj, rho_guess, box, mesh, components)
            box = initialize_box@dft_model(obj, rho_guess, box, mesh, components);
            
            box_bath = box.set_rho(box.rho0,mesh);
            
            mu_bath_excess = arrayfun(@(x)(x.calc_mu_bath(box_bath,mesh,components)),obj.excess,'UniformOutput',false);
            mu_bath_excess = cat(5,mu_bath_excess{:}); % creates an NX by NY by NZ by Ncomp by length(excess) array that stores the mu_bath's of the various excess functionals
                        
            mu_bath = obj.ideal_gas.calc_mu_bath(box_bath, mesh, components) + sum(mu_bath_excess,5);
            V_external = obj.external.calc_mu(box_bath, mesh, components);
            
            box = box.set_static_fields(mu_bath, V_external);
        end
        
        function [err, new_box] = calc_Picard_step(obj, box, mesh, components)
            mu_excess = arrayfun(@(x)(x.calc_mu(box,mesh,components)),obj.excess,'UniformOutput',false);
            mu_excess = cat(5,mu_excess{:}); % creates an NX by NY by NZ by Ncomp by length(excess) array that stores the mu's of the various excess functionals
            mu_excess = sum(mu_excess,5);
            
            new_rho = exp(box.mu_bath-(mu_excess + box.V_external));
            
            % This effectively makes the decision variable log(rho) instead
            % of rho in the Picard iteration. This is useful since log(rho)
            % is unbounded, but rho must be non-negative
            err = log(box.rho) - log(new_rho); 
            err = norm(err(:));

            if nargout > 1
                new_box = box.set_rho(new_rho,mesh);
            end
        end
        
    end
    
end

