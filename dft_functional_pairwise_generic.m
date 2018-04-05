classdef dft_functional_pairwise_generic < dft_functional_pairwise
    %dft_functional_pairwise_generic - 
    %   input a function handle and a parameter array and this will set up
    %   the k-space matrices in U{} needed for calculating mu and the free
    %   energy
    %
    % obj = dft_functional_pairwise_generic(u, parameters, mesh, components)
    % parameters : cell array of Ncomp by Ncomp matrices of input
    %   parameters, where each cell is a symmetric matrix of parameter
    %   values
    % u : function handle that accepts u(r,parameters_vec), to calculate 
    %   the potential function in real space, where parameters_vec is a
    %   vector the length of the parameters cell array
    %
    %     NOTE: this will also work for Oleksy & Hansen style correlation
    %     terms, but with u(r,parameters_vec) = -Delta_c(r,parameters_vec)
    %
    %     NOTE: all cutoffs and values inside the hard core is handled in u
    %     not in this class. This also means that there is no explicit 
    %     check that the potential will fit inside the defined mesh.



    methods
        function obj = dft_functional_pairwise_generic(u, parameters, mesh, components)
            % check type of u
            if ~isa(u,'function_handle')
                error('ERROR: u must be a function handle')
            end
            
            % check size and symmetry of parameters
            num_parameters = size(parameters);
            for p=1:num_parameters
                if ~issymmetric(parameters{p})
                    error('ERROR: parameters must be symmetric or composed of symmetric matrices')
                end
                if ~isequal(size(parameters{p}), [components.Ncomp, components.Ncomp])
                    error('ERROR: size of parameters does not match the number of types')
                end
            end
                        
            obj.U = cell(components.Ncomp);
            k_unique = unique(mesh.k);
            for a=1:components.Ncomp
                for b=a:components.Ncomp
                    dHS = (components.dHS(a)+components.dHS(b))/2;
                    
                    current_parameters = cellfun(@(x)x(a,b),parameters);
                    waypoints = unique([current_parameters, dHS]);
                    u_r = @(r)u(r,current_parameters);
                    U_unique = (4*pi./k_unique).*arrayfun(@(k)( integral(@(r)r.*u_r(r).*sin(k*r), 0, Inf, 'RelTol', 1e-12, 'Waypoints', waypoints ) ), k_unique);
                    U = zeros(size(mesh.k));
                    for n = 1:length(k_unique)
                        U(mesh.k==k_unique(n)) = U_unique(n);
                    end
                    
                    %the k=0 limit is calculated separately
                    U(1,1,1) = 4*pi*integral(@(r)r.^2.*u_r(r), 0, Inf, 'RelTol', 1e-12, 'Waypoints', waypoints );
                    
                    obj.U{a,b} = U;
                    obj.U{b,a} = U;
                end
            end
        end
        
        
        
    end
    
end


