classdef dft_functional_poisson < dft_functional
    %dft_functional_external
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        l_B = 0; % Bjerrum length, the distance at which the electrostatic energy between two unit charge ions is 1kT 
    end
    
    methods
        function obj = dft_functional_poisson(l_B, mesh, components)
            % FIXME
        end
        
        function F = calc_free_energy(obj, box, mesh, components)
            % FIXME
            F = sum(F(:));
        end
        
        function mu = calc_mu(obj, box, mesh, components)
            mu = %FIXME
        end
        
        % calculates the electrostatic potential by solving Poisson's
        % equation in k-space
        % Poisson's equation rearranged to use l_B and in units where
        % "charge" is the charge number
        % In the end the equation looks like:
        % del^2(beta * e * psi) = -4pi * l_B * sum_i( charge_i * rho_i )
        % where beta=1/kT, e=electron charge, psi=electrostatic potential
        function beta_e_psi = calc_psi(obj, box, mesh, components)
            charge_density = sum(reshape(components.charge,1,1,1,components.Ncomp).*box.rho,4);
            net_charge = sum(charge_density(:));
            if abs(net_charge) > 1e-14
                error('net charge is too large')
            end
            
            rhs = -4*pi*obj.l_B*chrage_density; % right hand side of Poisson's equation
            RHS = mesh.forward_fft(rhs); % rhs transformed to k-space
            
            % the solution to Poisson's equation in k-space is -RHS/k^2
            BETA_E_PSI = -RHS./(mesh.k.^2);
            % correct for the k=0 limit, this sets the integral of the
            % electrostaic potential, this can be set to anything
            % (hopefully)
            BETA_E_PSI(1) = 0; 
            
            % transform the solution back to real space
            beta_e_psi = mesh.inverse_fft(BETA_E_PSI);
        end
    end
    
end
