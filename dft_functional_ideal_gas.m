classdef dft_functional_ideal_gas < dft_functional
    %dft_functional_ideal_gas
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function F = calc_free_energy(obj, box, mesh, components)
            if isa(box,'dft_box_isaft')
                F = sum(box.weighted_rho_seg.*(log(box.rho_seg)-1),4); % FIXME: why do we ignore the weighting here?
            else
                F = sum(box.weighted_rho.*(log(box.rho)-1),4); % FIXME: why do we ignore the weighting here?
            end
            F = F.*mesh.cell_weight*mesh.DX^(mesh.D);
            F = sum(F(:));
        end
        
        function mu = calc_mu(obj, box, mesh, components)
            mu = log(box.rho); % FIXME: why do we ignore the weighting here?
        end

    end
    
end

