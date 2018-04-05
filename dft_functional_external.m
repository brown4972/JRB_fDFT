classdef dft_functional_external < dft_functional
    %dft_functional_external
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        v_ext = []; % external potential; NX by NY by NZ by Ncomp array
    end
    
    methods
        function obj = dft_functional_external(v_ext, mesh, components)
            if isequal(size(v_ext), [mesh.NX, mesh.NY, mesh.NZ, components.Ncomp]) || (isequal(size(v_ext), [mesh.NX, mesh.NY]) && components.Ncomp==1) || (isequal(size(v_ext), [mesh.NX, mesh.NY, mesh.NZ]) && components.Ncomp==1)
                obj.v_ext = v_ext;
            else
                error('ERROR: size of v_ext does not match the mesh and number of components')
            end
        end
        
        function F = calc_free_energy(obj, box, mesh, components)
            mu = obj.calc_mu(box, mesh, components);
            F = sum(box.weighted_rho.*mu,4); % weighted rho here to account for walls and BCs
            F = F.*mesh.cell_weight*mesh.DX^(mesh.D);
            F = sum(F(:));
        end
        
        function mu = calc_mu(obj, box, mesh, components)
            mu = obj.v_ext;
        end
    end
    
end
