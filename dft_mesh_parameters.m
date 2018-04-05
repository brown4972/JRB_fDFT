classdef dft_mesh_parameters
    % dft_box_properties - class that stores the parameters for the computational mesh for DFT
    %
    
    properties (SetAccess = immutable)
        D = 0; % dimensionality of space
        BCs = ''; % Boundary conditions; string, either 'p' or 'r' for periodic and reflective
        DX = 0; % mesh spacing (note should be a even fraction of the HS diameters)
        
        cell_weight = []; % NX by NY by NZ array of weights to account for BCs in various calculations
        
        NX = 1; % number of cells in X direction
        NY = 1; % number of cells in Y direction
        NZ = 1; % number of cells in Z direction
        
        x = []; % NX by NY by NZ array of x locations
        y = []; % NX by NY by NZ array of y locations
        z = []; % NX by NY by NZ array of z locations
        s = []; % NX by NY by NZ array of cylindrical radii from the z axis NOTE: this accounts for BCs
        r = []; % NX by NY by NZ array of spherical radii from the origin NOTE: this accounts for BCs
        
        NKX = 1; % number of cells in KX direction
        NKY = 1; % number of cells in KY direction
        NKZ = 1; % number of cells in KZ direction
        
        kx = []; % NKX by NKY by NKZ array of x components of the k vector
        ky = []; % NKX by NKY by NKZ array of y components of the k vector
        kz = []; % NKX by NKY by NKZ array of z components of the k vector
        k = []; % NKX by NKY by NKZ array of k vector magnitudes
        
        forward_fft = @do_nothing; % function handle for forward FFT
        forward_fft_vx = @do_nothing; % function handle for forward FFT for the x component of the vector valued weighted densities (needed to account for symmetry in reflective BCs)
        forward_fft_vy = @do_nothing; % function handle for forward FFT for the y component of the vector valued weighted densities (needed to account for symmetry in reflective BCs)
        forward_fft_vz = @do_nothing; % function handle for forward FFT for the z component of the vector valued weighted densities (needed to account for symmetry in reflective BCs)
        inverse_fft = @do_nothing; % function handle for inverse FFT
    end
    
    methods
        function obj = dft_mesh_parameters(D, BCs, DX, NX, NY, NZ)
            obj.D = D;
            
            if length(DX) ~= 1
                error('ERROR: mesh spacing must be a single number')
            end
            obj.DX = DX;
            
            if length(NX) ~= 1
                error('ERROR: box size must be a single number')
            end
            obj.NX = NX;

            if D>=2
                if length(NY) ~= 1
                    error('ERROR: box size must be a single number')
                end
                obj.NY = NY;
            end
            
            if D==3
                if length(NZ) ~= 1
                    error('ERROR: box size must be a single number')
                end
                obj.NZ = NZ;
            end
            
            if BCs ~= 'p' && BCs ~= 'r'
                error('ERROR: BCs must be "p" or "r"')
            end
            obj.BCs = BCs;
                 
            obj.cell_weight = ones(obj.NX, obj.NY, obj.NZ);
            
            if BCs == 'p'
                obj.NKX = obj.NX;
                obj.NKY = obj.NY;
                obj.NKZ = obj.NZ;
                
                % periodic boundary conditions are automatically satisfied
                % by the assumtions of the standard discrete Fourier
                % transform
                switch D
                    case 1
                        obj.forward_fft = @fft;
                        obj.inverse_fft = @(x)ifft(x,'symmetric');
                    case 2
                        obj.forward_fft = @fft2;
                        obj.inverse_fft = @(x)ifft2(x,'symmetric');
                    case 3
                        obj.forward_fft = @fft3;
                        obj.inverse_fft = @(x)ifft3(x,'symmetric');
                end
                
                obj.forward_fft_vx = obj.forward_fft;
                obj.forward_fft_vy = obj.forward_fft;
                obj.forward_fft_vz = obj.forward_fft;
                
            else %if BCs == 'r'
                switch D
                    case 1
                        obj.NKX = 2*obj.NX-2;
                        
                        obj.cell_weight(1,1,1) = 0.5;
                        obj.cell_weight(obj.NX,1,1) = 0.5;
                        
                        % The definition of the forward FFT duplicates the
                        % data to appropriately apply reflective symmetry,
                        % NOTE that the vector data is both flipped left 
                        % to right and made negative 
                        % NOTE also that the k-space arrays are now twice
                        % as large as the real-space ones
                        obj.forward_fft = @(x)real(fft([x;flip(x(2:obj.NX-1,1,1,:),1)]));
                        obj.forward_fft_vx = @(x)fft([x;-flip(x(2:obj.NX-1,1,1,:),1)]);
                        obj.forward_fft_vy = obj.forward_fft;
                        obj.forward_fft_vz = obj.forward_fft;
                        
                        % This one-liner throws out the extra symmetrical
                        % data that was added in the forward FFT
                        obj.inverse_fft = @(X)subsref(ifft(X,'symmetric'),struct('type','()', 'subs', {{1:obj.NX,1,1,':'}}));

                    case 2
                        warning('2D reflective BCs are buggy')
                        
                        obj.NKX = 2*obj.NX-2;
                        obj.NKY = 2*obj.NY-2;
                        
                        obj.cell_weight(1,:,1) = 0.5;
                        obj.cell_weight(:,1,1) = 0.5;
                        obj.cell_weight(obj.NX,:,1) = 0.5;
                        obj.cell_weight(:,obj.NY,1) = 0.5;
                        obj.cell_weight(1,1,1) = 0.25;
                        obj.cell_weight(obj.NX,1,1) = 0.25;
                        obj.cell_weight(1,obj.NY,1) = 0.25;
                        obj.cell_weight(obj.NX,obj.NY,1) = 0.25;
                        
                        % The definition of the forward FFT duplicates the
                        % data to appropriately apply reflective symmetry,
                        % NOTE that the vector data is both flipped left to
                        % right or up to down and made negative as 
                        % appripriate for the particular quadrant 
                        % NOTE also that the k-space arrays are now 4 times
                        % as large as the real-space ones
                        obj.forward_fft = @(x)real(fft2([x, flip(x(1:obj.NX,2:obj.NY-1,1,:),2); flip(x(2:obj.NX-1,1:obj.NY,1,:),1), flip(flip(x(2:obj.NX-1,2:obj.NY-1,1,:),1),2)]));
                        obj.forward_fft_vx = @(x)fft2([x, flip(x(1:obj.NX,2:obj.NY-1,1,:),2); -flip(x(2:obj.NX-1,1:obj.NY,1,:),1), -flip(flip(x(2:obj.NX-1,2:obj.NY-1,1,:),1),2)]);
                        obj.forward_fft_vy = @(x)fft2([x, -flip(x(1:obj.NX,2:obj.NY-1,1,:),2); flip(x(2:obj.NX-1,1:obj.NY,1,:),1), -flip(flip(x(2:obj.NX-1,2:obj.NY-1,1,:),1),2)]);
                        obj.forward_fft_vz = obj.forward_fft;
                        
                        % This one-liner throws out the extra symmetrical
                        % data that was added in the forward FFT
                        obj.inverse_fft = @(X)subsref(ifft2(X,'symmetric'),struct('type','()', 'subs', {{1:obj.NX,1:obj.NY,1,':'}}));

                    case 3
                        error('3D reflective BCs not implemented')
                end
            end
            
            % set up references for real space coordinates
            obj.x = 0:obj.DX:(obj.NX-1)*obj.DX;
            obj.y = 0:obj.DX:(obj.NY-1)*obj.DX;
            obj.z = 0:obj.DX:(obj.NZ-1)*obj.DX;
            [obj.y,obj.x,obj.z] = meshgrid(obj.y,obj.x,obj.z); % NOTE: x and y are "reversed" here so that the x direction is the first index, as is done elsewhere 
            
            if BCs == 'p'
                obj.r = sqrt(min(min(obj.x.^2,(obj.x-obj.NX*obj.DX).^2),(obj.x+obj.NX*obj.DX).^2) ...
                            +min(min(obj.y.^2,(obj.y-obj.NY*obj.DX).^2),(obj.y+obj.NY*obj.DX).^2) ...
                            +min(min(obj.z.^2,(obj.z-obj.NZ*obj.DX).^2),(obj.z+obj.NZ*obj.DX).^2));
                obj.s = sqrt(min(min(obj.x.^2,(obj.x-obj.NX*obj.DX).^2),(obj.x+obj.NX*obj.DX).^2) ...
                            +min(min(obj.y.^2,(obj.y-obj.NY*obj.DX).^2),(obj.y+obj.NY*obj.DX).^2));
            else %if BCs == 'r'
                obj.r = sqrt(obj.x.^2 + obj.y.^2 + obj.z.^2);
                obj.s = sqrt(obj.x.^2 + obj.y.^2);
            end
            
            % set up references for k-space coordinates
            obj.kx = 0:obj.NKX-1;
            obj.ky = 0:obj.NKY-1;
            obj.kz = 0:obj.NKZ-1;
            
            obj.kx(0:obj.NKX-1>obj.NKX/2) = obj.kx(0:obj.NKX-1>obj.NKX/2) - obj.NKX;
            obj.ky(0:obj.NKY-1>obj.NKY/2) = obj.ky(0:obj.NKY-1>obj.NKY/2) - obj.NKY;
            obj.kz(0:obj.NKZ-1>obj.NKZ/2) = obj.kz(0:obj.NKZ-1>obj.NKZ/2) - obj.NKZ;
            
            obj.kx = obj.kx*2*pi/obj.DX/obj.NKX;
            obj.ky = obj.ky*2*pi/obj.DX/obj.NKY;
            obj.kz = obj.kz*2*pi/obj.DX/obj.NKZ;
            
            [obj.ky, obj.kx, obj.kz] = meshgrid(obj.ky,obj.kx,obj.kz); % NOTE: same x/y reversal as above
            
            obj.k = sqrt(obj.kx.^2 + obj.ky.^2 + obj.kz.^2);
        end
    end
    
end

