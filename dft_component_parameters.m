classdef dft_component_parameters
    % dft_component_properties - class that stores the parameters of components for a DFT system
    %
    
    properties (SetAccess = protected)
        Ncomp = 0; % number of component types in the system
        Nseg = 0; % total number of segments in the system (same as Ncomp for non-polymer systems)
        rho0 = []; % bath density; length Ncomp array
        dHS = []; % particle hard sphere diameters; length Ncomp array
        seg_to_type = []; % vector that maps segment number to its type; length Nseg vector
    end
    
    methods
        function obj = dft_component_parameters(Ncomp, dHS, rho0)
            if length(rho0) == Ncomp
                obj.Ncomp = Ncomp;
                obj.rho0 = rho0;
            else
                error('ERROR: Bath densities do not match number of components')
            end
                        
            if length(dHS) == Ncomp
                obj.dHS = dHS;
            else
                error('ERROR: Hard sphere diameters do not match number of components')
            end
            
            obj.Nseg = Ncomp;
            obj.seg_to_type = 1:Ncomp;
            
        end
    end
    
end

