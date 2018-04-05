classdef dft_model
    % dft_model - base class to store functionals and calculate Omega and
    %   NL solver update steps
    %   This is the base model class - a subclass is needed to do specifics
    
    properties (SetAccess = protected)
        ideal_gas = dft_functional; % ideal gas functional
        external = dft_functional; % external field functional
        excess = dft_functional; % vector of excess functionals
    end
    
    methods
        function obj = dft_model(ideal_gas, external, excess, mesh, components)
            if ~isa(ideal_gas,'dft_functional')
                error('ERROR: ideal_gas must be a dft_functional')
            elseif ~isa(external,'dft_functional')
                error('ERROR: external must be a dft_functional')
            elseif ~isa(excess,'dft_functional')
                error('ERROR: external must be composed of type dft_functional')
            elseif ~isvector(excess)
                error('ERROR: external must be a vector')
            end
            
            obj.ideal_gas = ideal_gas;
            obj.external = external;
            obj.excess = excess;
        end
        
        function box = initialize_box(obj, rho_guess, box, mesh, components)
            box = box.set_rho(rho_guess, mesh, components);
        end
        
        function [Omega, Omega_ideal_gas, Omega_external, Omega_mu, Omega_excess] = calc_Omega(obj, box, mesh, components)            
            Omega_ideal_gas = obj.ideal_gas.calc_free_energy(box, mesh, components);
            Omega_external = obj.external.calc_free_energy(box, mesh, components);
            Omega_mu = calc_Omega_mu(obj, box, mesh, components);
            Omega_excess = arrayfun(@(x)(x.calc_free_energy(box,mesh,components)),obj.excess);
            Omega = Omega_ideal_gas + Omega_external + Omega_mu + sum(Omega_excess);
        end
        
        % Omega_mu is the contribution to the grand free energy from the
        % fixed bath chemical potential
        function [ Omega_mu ] = calc_Omega_mu(obj, box, mesh, components)
            Omega_mu = -sum(box.weighted_rho.*box.mu_bath,4);
            Omega_mu = Omega_mu.*mesh.cell_weight*mesh.DX^(mesh.D);
            Omega_mu = sum(Omega_mu(:));
        end
        
        function [err,new_box] = calc_Picard_step(obj, box, mesh, components)
            err = 0;
            
            if nargout > 1
                new_box = box;
            end
        end
        
        function new_box = calc_Newton_step(obj, box, mesh, components)
            [~, new_box] = obj.calc_Picard_step(box, mesh, components);
        end
    end
    
end

