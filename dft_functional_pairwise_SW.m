classdef dft_functional_pairwise_SW < dft_functional_pairwise
    %dft_functional_pairwise_SW - square well potential
    %   uses closed form solutions to the transforms of the SW potential
    %   into k-space
    
    properties (SetAccess = immutable)
        epsilon = []; % Ncomp by Ncomp symmetric matrix of epsilon_ij values (note a negative epsilon means a repulsion)
        sigma = []; % Ncomp by Ncomp symmetric matrix of sigma_ij values
        r_cut = []; % Ncomp by Ncomp symmetric matrix of r_cut_ij values
    end
    
    methods
        function obj = dft_functional_pairwise_SW(epsilon, r_cut, mesh, components)
            if isequal(size(epsilon),[components.Ncomp, components.Ncomp]) && issymmetric(epsilon)
                obj.epsilon = epsilon;
            else
                error('ERROR: size of epsilon does not match the number of types or is not symmetric')
            end
            
            if isequal(size(r_cut),[components.Ncomp, components.Ncomp]) && issymmetric(r_cut)
                obj.r_cut = r_cut;
            else
                error('ERROR: size of r_cut does not match the number of types or is not symmetric')
            end
            
            r_c_max = max(r_cut(:));
            N_XYZ = [mesh.NX, mesh.NY, mesh.NZ];
            L_min = mesh.DX*min(N_XYZ(N_XYZ>1));
            if 2*r_c_max > L_min
                error('ERROR: 2*max(r_cut) is larger than the smallest box dimension') % FIXME: this could be accounted for, would need to change the way u(r) is calculated
            end
            
            
            obj.U = cell(components.Ncomp);
            for a=1:components.Ncomp
                for b=a:components.Ncomp
                    eps = obj.epsilon(a,b);
                    if eps ~= 0
                        dHS = (components.dHS(a)+components.dHS(b))/2;
                        r_c = obj.r_cut(a,b);
                        obj.U{a,b} = -eps*(4*pi./mesh.k.^3).*( (sin(r_c.*mesh.k)-r_c.*mesh.k.*cos(r_c.*mesh.k)) - (sin(dHS.*mesh.k)-dHS.*mesh.k.*cos(dHS.*mesh.k)) );
                        obj.U{a,b}(1,1,1) = -eps*(4/3)*pi*(r_c.^3 - dHS.^3); % correct for k=0 limit
                    else
                        obj.U{a,b} = 0;
                    end
                    if a~=b
                        obj.U{b,a} = obj.U{a,b};
                    end
                end
            end
        end
        
    end
    
end


