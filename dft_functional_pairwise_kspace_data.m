classdef dft_functional_pairwise_kspace_data < dft_functional_pairwise
    %dft_functional_pairwise_generic - 
    %   input k-space data for a pairwise interaction and this will set up
    %   the k-space matrices in U{} needed for calculating mu and the free
    %   energy, interpolating with splines, where necessary
    %
    % obj = dft_functional_pairwise_kspace_data(k_in, U_in, mesh, components)
    % k_in: vector is k-values that correspond to the U_in values
    %
    % U_in: Ncomp by Ncomp cell array of vectors representing the k-space
    %   values of a pairwise interaction
    %
    %     NOTE: this will also work for Oleksy & Hansen style correlation
    %     terms, but with U_in = -Delta_C(k)
    %
    %     NOTE: if max(k_in) < max(mesh.k), all U{a,b} values at k >
    %     max(k_in) are set to zero. min(k_in) < min(mesh.k) == 0, the
    %     U{a,b} values are extrapolated.



    methods
        function obj = dft_functional_pairwise_kspace_data(k_in, U_in, mesh, components)
            % check size of U_in
            if ~isequal(size(U_in), [components.Ncomp, components.Ncomp])
                error('ERROR: size of U_in does not match the number of types')
            end
                
            % check symmetry of U_in, and that each entry matches k_in
            for a=1:components.Ncomp
                for b=a:components.Ncomp
                    if ~isequal(size(k_in), size(U_in{a,b}))
                        error('ERROR: size of every entry of U_in does not match the size of k_in')
                    end
                    if ~isequal(U_in{a,b}, U_in{b,a})
                        error('ERROR: U_in must be symmetric')
                    end
                end
            end
                        
            obj.U = cell(components.Ncomp);
            k_unique = unique(mesh.k);
            for a=1:components.Ncomp
                for b=a:components.Ncomp
%                     dHS = (components.dHS(a)+components.dHS(b))/2;
                    
                    U_unique = interp1(k_in,U_in{a,b},k_unique,'spline','extrap');
                    U_unique(k_unique>max(k_in)) = 0;
                    
                    U = zeros(size(mesh.k));
                    for n = 1:length(k_unique)
                        U(mesh.k==k_unique(n)) = U_unique(n);
                    end
                    
                    obj.U{a,b} = U;
                    obj.U{b,a} = U;
                end
            end
        end
        
        
        
    end
    
end


