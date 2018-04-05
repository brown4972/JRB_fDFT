classdef dft_functional < matlab.mixin.Heterogeneous
    %dft_functional
    %   Base class for DFT functionals, returns a zero free energy and
    %   (inhomogeneous) mu if used by itself, also handles the homogeneous 
    %   "bath" mu for functionals where it is not overridden
    
    properties
    end
    
    methods
        function F = calc_free_energy(obj, box, mesh, components)
            F = 0;
        end
        
        function mu = calc_mu(obj, box, mesh, components)
            mu = zeros(size(box.rho));
        end
        
        function mu_bath = calc_mu_bath(obj, box, mesh, components)
            box0 = box.set_rho(box.rho0,mesh);
            mu_bath = obj.calc_mu(box0, mesh, components);
        end
        
    end
    
end

