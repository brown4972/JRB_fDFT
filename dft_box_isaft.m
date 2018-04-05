classdef dft_box_isaft < dft_box
    %dft_box_isaft Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        rho_seg = []; % segmental densities, NX by NX by NZ by Nseg array
        weighted_rho_seg = []; % segmental densities weighted to account for walls, NX by NX by NZ by Nseg array
        
        % nonlocal densities for iSAFT
        xi2 = []; % NX by NX by NZ array
        xi3 = []; % NX by NX by NZ array
    end
    
    properties (SetAccess = immutable)
        % weighting functions for iSAFT nonlocal densities, in k-space:
        WXI2 = []; % NX by NX by NZ by Ncomp array
        WXI3 = []; % NX by NX by NZ by Ncomp array
        
        % weighting functions for iSAFT bonds, in k-space:
        Wbond = []; % NX by NX by NZ by Nseg array
        
        % weighting functions for iSAFT associations (breakable bonds), in k-space:
        Wassoc = []; % NX by NX by NZ by Ncomp by Ncomp array
    end
    
    
    methods
        function obj = dft_box_isaft(mesh, components, wall_edge_location)
            if ~isa(components,'dft_polymer_parameters')
                error('ERROR: components must be a dft_polymer_parameters')
            end
            
            if nargin <= 2
                wall_edge_location = logical(zeros(mesh.NX,mesh.NY,mesh.NZ));
            end
            obj = obj@dft_box(mesh, components, wall_edge_location);

            
            obj.mu_bath = zeros(mesh.NX, mesh.NY, mesh.NZ,components.Npoly);
            
            obj.rho_seg = zeros(mesh.NX, mesh.NY, mesh.NZ,components.Nseg);
            obj.weighted_rho_seg = zeros(mesh.NX, mesh.NY, mesh.NZ,components.Nseg);
            
            % weighting functions for the iSAFT nonlocal densities
            dHS = reshape(components.dHS,1,1,1,components.Ncomp);
            obj.WXI2 = (pi./mesh.k.^3).*(sin(dHS.*mesh.k)-dHS.*mesh.k.*cos(dHS.*mesh.k))./(2*dHS);
            obj.WXI2(1,1,1,:) = (1/6)*pi*dHS.^2; % correct for k=0 limit
            obj.WXI3 = obj.WXI2.*dHS;
            
            
            % weighting functions for iSAFT bonds
            dBond = zeros(1,1,1,components.Nseg);
            for seg=2:components.Nseg
                if components.seg_to_poly(seg) == components.seg_to_poly(seg-1)
                    d_seg = components.dHS(components.seg_to_type(seg));
                    d_prev = components.dHS(components.seg_to_type(seg-1));
                    dBond(seg) = (d_seg+d_prev)/2;
                end
            end
            obj.Wbond = sin(dBond.*mesh.k)./(dBond.*mesh.k);
            obj.Wbond(1,1,1,:) = 1; % correct for k=0 limit
            obj.Wbond(:,:,:,dBond==0) = 0; % set the nonbonded terms weigting functions to zero
            
            
            % weighting functions for iSAFT associations (breakable bonds)
            dAssoc = zeros(1,1,1,components.Ncomp,components.Ncomp);
            for m=1:components.Ncomp
                for n=1:components.Ncomp
                    dAssoc(1,1,1,m,n) = (components.dHS(m)+components.dHS(n))/2;
                end
            end
            obj.Wassoc = (4*pi*dAssoc).*sin(dAssoc.*mesh.k)./(mesh.k);
            obj.Wassoc(1,1,1,:,:) = 4*pi*dAssoc.^2; % correct for k=0 limit
        end
        
        % can set densities either from explicit segmental densities, or by
        % component densities (which are then scaled to make segmental
        % densities)
        function obj = set_rho(obj, new_rho, mesh, components)
            if isequal( size(new_rho), size(obj.rho_seg))
                obj.rho_seg = new_rho;
                new_rho = zeros(size(obj.rho));
                for n=1:components.Ncomp
                    new_rho(:,:,:,n) = sum(obj.rho_seg(:,:,:,components.seg_to_type == n),4);
                end
            elseif isequal( size(new_rho), size(obj.rho) )
                for seg=1:components.Nseg
                    obj.rho_seg(:,:,:,seg) = new_rho(:,:,:,components.seg_to_type(seg))/sum(components.seg_to_type == components.seg_to_type(seg));
                end
            else
                error('ERROR: the size of new_rho must match rho or rho_seg')
            end
                 
            % only apply wall weighting when the density isn't bulk (this is a bit of a hack because not all functionals have separate mu_bulk calculations but instead just use the inhomogeneous one)
            if isequal(obj.rho,obj.rho0)
                obj.weighted_rho_seg = obj.rho_seg;
            else
                obj.weighted_rho_seg = obj.rho_seg.*obj.wall_weight(:,:,:,components.seg_to_type);
            end
            
            obj = set_rho@dft_box(obj, new_rho, mesh, components);
            
            XI2 = sum(obj.RHO.*obj.WXI2,4);
            XI3 = sum(obj.RHO.*obj.WXI3,4);
            
            obj.xi2 = mesh.inverse_fft(XI2);
            obj.xi3 = mesh.inverse_fft(XI3);
        end
        
        % the mix routine mixes the densities and nonlocal densities (as
        % well as the fields) of two boxes, to avoid having to recompute
        % the nonlocal densities
        function obj = mix(obj, obj2, alpha)
            if ~isa(obj2,'dft_box_isaft')
                error('ERROR: must input objects of type dft_box_isaft')
            end
            
            obj = mix@dft_box(obj, obj2, alpha);
            
            obj.rho_seg = (1-alpha)*obj.rho_seg + alpha*obj2.rho_seg;
            
            obj.xi2 = (1-alpha)*obj.xi2 + alpha*obj2.xi2;
            obj.xi3 = (1-alpha)*obj.xi3 + alpha*obj2.xi3;
        end

        % helper functions for converting back and forth into the decision
        % variable (log(rho)) vector used in the Newton solver
        function v = vec(obj)
            v = log(obj.rho_seg(:));
        end
        
        function obj = unvec(obj,v,mesh,components)
            obj = obj.set_rho(reshape(exp(v),mesh.NX,mesh.NY,mesh.NZ,components.Nseg), mesh, components);
        end
    end
    
end

