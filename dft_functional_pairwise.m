classdef dft_functional_pairwise < dft_functional
    %dft_functional_pairwise - base class for pairwise interactions
    %   holds the generic calc_mu and calc_free_energy routines, but they
    %   will throw an error withou a subclass to define the U{} matrices. 
    
    properties (SetAccess = protected)
        U = {}; % cell array of the k-space potential functions
    end
    
    methods
        
        function F = calc_free_energy(obj, box, mesh, components)
            mu = obj.calc_mu(box, mesh, components);
            F = 0.5*sum(box.weighted_rho.*mu,4);
            F = F.*mesh.cell_weight*mesh.DX^(mesh.D);
            F = sum(F(:));
        end
        
        function mu = calc_mu(obj, box, mesh, components)
            % FIXME: vectorize?
            mu = zeros(mesh.NKX, mesh.NKY, mesh.NKZ, components.Ncomp);
            for a=1:components.Ncomp
                for b=1:components.Ncomp
                    mu(:,:,:,a) = mu(:,:,:,a) + obj.U{a,b}.*box.RHO(:,:,:,b);
                end
            end
            
            mu = mesh.inverse_fft(mu);
        end
        
    end
    
end


