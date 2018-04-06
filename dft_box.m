classdef dft_box
    %dft_box Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        rho = []; % density profiles; NX by NX by NZ by Ncomp array
        weighted_rho = []; % density profiles weighted to account for walls; NX by NX by NZ by Ncomp array
        RHO = []; % WEIGHTED density profiles in k-space; NX by NX by NZ by Ncomp array
                
        % weighted densities for hard sphere functionals:
        n0 = []; % NX by NX by NZ array
        n1 = []; % NX by NX by NZ array
        n2 = []; % NX by NX by NZ array
        n3 = []; % NX by NX by NZ array
        nv1x = []; % NX by NX by NZ array
        nv1y = []; % NX by NX by NZ array
        nv1z = []; % NX by NX by NZ array
        nv2x = []; % NX by NX by NZ array
        nv2y = []; % NX by NX by NZ array
        nv2z = []; % NX by NX by NZ array
        
        % static fields to carry around
        rho0 = []; %bath density; NX by NX by NZ by Ncomp array
        mu_bath = []; % bath chemical potential; NX by NX by NZ by Ncomp (or Npoly for polymers) array
        V_external = []; % external field; NX by NX by NZ by Ncomp array
        wall_weight = []; % array to weight densities when they are adjacent to a wall; NX by NX by NZ by Ncomp array
    end
        
    properties (SetAccess = immutable)
        % weighting functions for hard sphere weighted densities, in k-space:
        W0 = []; % NX by NX by NZ by Ncomp array
        W1 = []; % NX by NX by NZ by Ncomp array
        W2 = []; % NX by NX by NZ by Ncomp array
        W3 = []; % NX by NX by NZ by Ncomp array
        Wv1x = []; % NX by NX by NZ by Ncomp array
        Wv1y = []; % NX by NX by NZ by Ncomp array
        Wv1z = []; % NX by NX by NZ by Ncomp array
        Wv2x = []; % NX by NX by NZ by Ncomp array
        Wv2y = []; % NX by NX by NZ by Ncomp array
        Wv2z = []; % NX by NX by NZ by Ncomp array
    end
    
    
    methods
        function obj = dft_box(mesh, components, wall_edge_location)
            if ~isa(mesh,'dft_mesh_parameters') || ~isa(components, 'dft_component_parameters')
                error('ERROR: must input objects of dft_mesh_parameters and dft_component_parameters')
            end
            
            obj.rho0 = zeros(mesh.NX, mesh.NY, mesh.NZ,components.Ncomp);
            for n=1:components.Ncomp
                obj.rho0(:,:,:,n) = components.rho0(n);
            end
            obj.rho = obj.rho0;
            
            obj.mu_bath = zeros(mesh.NX, mesh.NY, mesh.NZ,components.Ncomp);
            
            % set up half-cell weights for walls -- to get walls to be
            % "halfway" though their cell, these cells are wighted by half.
            % Corners are weighted by 1/4 or 1/8 depending on geometry
            obj.wall_weight = ones(1,1,1,components.Ncomp);
            if nargin >= 3 && sum(sum(sum(wall_edge_location))) > 0
                obj.wall_weight = ones(mesh.NX, mesh.NY, mesh.NZ,components.Ncomp);
                if ~islogical(wall_edge_location)
                    error('wall_edge_location must be a logical array')
                end
                
                if ~isequal(size(wall_edge_location), size(mesh.cell_weight))
                    error('wall_edge_location must have the same dimensions as the box')
                end
                
                switch mesh.D
                    case 1
                        for n=1:components.Ncomp
                            obj.wall_weight(circshift(wall_edge_location, components.dHS(n)/mesh.DX/2),n) = 0.5;
                            obj.wall_weight(circshift(wall_edge_location,-components.dHS(n)/mesh.DX/2),n) = 0.5;
                        end
                    case 2                    
                        pwr =  circshift(wall_edge_location,  components.dHS(n)/mesh.DX/2,1) + ... 
                               circshift(wall_edge_location, -components.dHS(n)/mesh.DX/2,1) + ...
                               circshift(wall_edge_location,  components.dHS(n)/mesh.DX/2,2) +  ...
                               circshift(wall_edge_location, -components.dHS(n)/mesh.DX/2,2);
                        obj.wall_weight = 2.^(-pwr);
                    case 3
                        error('3D walls not yet implemented')
                end
            end
            
            % weighting functions for the hard sphere weighted densities
            R = reshape(components.dHS/2,1,1,1,components.Ncomp);
            obj.W0 = sin(R.*mesh.k)./(R.*mesh.k);
            obj.W0(1,1,1,:) = 1; % correct for k=0 limit
            obj.W1 = obj.W0.*R;
            obj.W2 = obj.W1.*(4*pi*R);
            obj.W3 = (4*pi./mesh.k.^3).*(sin(R.*mesh.k)-R.*mesh.k.*cos(R.*mesh.k));
            obj.W3(1,1,1,:) = (4/3)*pi*R.^3; % correct for k=0 limit
            
            obj.Wv2x = -1i*obj.W3.*mesh.kx;
            obj.Wv2y = -1i*obj.W3.*mesh.ky;
            obj.Wv2z = -1i*obj.W3.*mesh.kz;
            
            obj.Wv1x = obj.Wv2x./(4*pi*R);
            obj.Wv1y = obj.Wv2y./(4*pi*R);
            obj.Wv1z = obj.Wv2z./(4*pi*R);
        end
        
        function obj = set_rho(obj, new_rho, mesh, components)
            if ~isequal(size(new_rho), size(obj.rho))
                error('new_rho must be the same size as the current rho')
            end
            
            obj.rho = new_rho;
            
            % only apply weighting when the density isn't bulk
            % FIXME: this is a bit of a hack because not all functionals have separate mu_bulk calculations but instead just use the inhomogeneous one
            if isequal(obj.rho,obj.rho0)
                obj.weighted_rho = obj.rho;
            else
                obj.weighted_rho = obj.rho.*obj.wall_weight; % apply wall weighting here
            end
            
            obj.RHO = mesh.forward_fft(obj.weighted_rho); % using the wighted density here automatically applies the weighting for the HS and pairwise functionals
            
            % Calculate the weighted densities
            % first calculate in k-space
            N0 = sum(obj.RHO.*obj.W0,4);
            N1 = sum(obj.RHO.*obj.W1,4);
            N2 = sum(obj.RHO.*obj.W2,4);
            N3 = sum(obj.RHO.*obj.W3,4);
            
            Nv1x = sum(obj.RHO.*obj.Wv1x,4);
            Nv1y = sum(obj.RHO.*obj.Wv1y,4);
            Nv1z = sum(obj.RHO.*obj.Wv1z,4);
            
            Nv2x = sum(obj.RHO.*obj.Wv2x,4);
            Nv2y = sum(obj.RHO.*obj.Wv2y,4);
            Nv2z = sum(obj.RHO.*obj.Wv2z,4);
            
            % inverse transforms
            obj.n0 = mesh.inverse_fft(N0);
            obj.n1 = mesh.inverse_fft(N1);
            obj.n2 = mesh.inverse_fft(N2);
            obj.n3 = mesh.inverse_fft(N3);
            
            obj.nv1x = mesh.inverse_fft(Nv1x);
            obj.nv1y = mesh.inverse_fft(Nv1y);
            obj.nv1z = mesh.inverse_fft(Nv1z);
            
            obj.nv2x = mesh.inverse_fft(Nv2x);
            obj.nv2y = mesh.inverse_fft(Nv2y);
            obj.nv2z = mesh.inverse_fft(Nv2z);
        end
        
        function obj = set_static_fields(obj, mu_bath, V_external)
            if ~isequal(size(V_external), size(obj.rho))
                error('V_external must be the same size as the current rho')
            elseif ~isequal(size(mu_bath), size(obj.mu_bath))
                error('mu_bath must be NX by NX by NZ by Ncomp (or Npoly for polymers)')
            end
            
            obj.mu_bath = mu_bath;
            obj.V_external = V_external;
        end
        
        % the mix routine mixes the densities and nonlocal densities of two
        % boxes, to avoid having to recompute the nonlocal densities
        function obj = mix(obj, obj2, alpha)
            if ~isa(obj2,'dft_box')
                error('ERROR: must input objects of type dft_box')
            elseif 0 >= alpha || alpha >= 1
                error('ERROR: alpha must be between 0 and 1')
            end
            
            obj.rho = (1-alpha)*obj.rho + alpha*obj2.rho;
            obj.weighted_rho = (1-alpha)*obj.weighted_rho + alpha*obj2.weighted_rho;
            obj.RHO = (1-alpha)*obj.RHO + alpha*obj2.RHO;
            
            obj.n0 = (1-alpha)*obj.n0 + alpha*obj2.n0;
            obj.n1 = (1-alpha)*obj.n1 + alpha*obj2.n1;
            obj.n2 = (1-alpha)*obj.n2 + alpha*obj2.n2;
            obj.n3 = (1-alpha)*obj.n3 + alpha*obj2.n3;
            
            obj.nv1x = (1-alpha)*obj.nv1x + alpha*obj2.nv1x;
            obj.nv1y = (1-alpha)*obj.nv1y + alpha*obj2.nv1y;
            obj.nv1z = (1-alpha)*obj.nv1z + alpha*obj2.nv1z;
            
            obj.nv2x = (1-alpha)*obj.nv2x + alpha*obj2.nv2x;
            obj.nv2y = (1-alpha)*obj.nv2y + alpha*obj2.nv2y;
            obj.nv2z = (1-alpha)*obj.nv2z + alpha*obj2.nv2z;
        end
        
        % helper functions for converting back and forth into the decision
        % variable (log(rho)) vector used in the Newton solver
        function v = vec(obj)
            v = log(obj.rho(:));
        end
        
        function obj = unvec(obj,v,mesh,components)
            obj = set_rho(obj, reshape(exp(v),mesh.NX,mesh.NY,mesh.NZ,components.Ncomp), mesh);
        end
    end
end

