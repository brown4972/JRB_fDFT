classdef dft_polymer_parameters < dft_component_parameters
    % dft_component_properties - class that stores the polymer parameters
    %   for linear polymers only
    
    properties (SetAccess = protected)
        Npoly = 0; % number of polymer types in the system
        rho_poly_0 = []; % bath density of polymers; length Npoly array
        poly_len = []; % length of polymers; length Npoly vector
        sequences = {} % cell array of polymer sequences (i.e. their component types); length Npoly cell array of vectors
        nGamma = []; % number of bonds per segment; length Nseg vector
        seg_to_poly = []; % mapping from segment number to polymer number; length Nseg vector
    end
    
    methods
        function obj = dft_polymer_parameters(Ncomp, dHS, Npoly, rho_poly_0, sequences)
            obj = obj@dft_component_parameters(Ncomp, dHS, zeros(1,Ncomp));
            
            if length(sequences) ~= Npoly
                error('ERROR: sequences must have length Npoly')
            elseif length(rho_poly_0) ~= Npoly
                error('ERROR: rho_poly_0 must have length Npoly')
            end
            
            obj.Npoly = Npoly;
            obj.rho_poly_0 = rho_poly_0;
            obj.sequences = sequences;
            
            obj.Nseg = sum(cellfun(@length,sequences));
            obj.rho0 = zeros(1,obj.Ncomp);
            obj.poly_len = zeros(1,Npoly);
            for poly=1:Npoly
                obj.poly_len(poly) = length(sequences{poly});
                for a=1:Ncomp
                    obj.rho0(a) = obj.rho0(a) + rho_poly_0(poly)*sum(sequences{poly}==a)/obj.poly_len(poly);
                end
            end
            
            obj.nGamma = zeros(1,obj.Nseg);
            obj.seg_to_poly = zeros(1,obj.Nseg);
            obj.seg_to_type = zeros(1,obj.Nseg);
            seg = 0;
            for poly = 1:Npoly
                for n = 1:obj.poly_len(poly)
                    seg = seg+1;
                    obj.nGamma(seg) = sum([n>1,n<obj.poly_len(poly)]);
                    obj.seg_to_poly(seg) = poly;
                    obj.seg_to_type(seg) = sequences{poly}(n);
                end
            end
                
            
        end
    end
    
end
