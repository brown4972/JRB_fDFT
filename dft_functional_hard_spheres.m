classdef dft_functional_hard_spheres < dft_functional
    % dft_functional_hard_spheres - base class for hard sphere functionals
    %   holds the generic calc_mu routine, but it will throw an error if
    %   the dPhi_dnX's are not defined (by a subclass)  
    
    properties (Access = protected)
        dPhi_dn0 = [];
        dPhi_dn1 = [];
        dPhi_dn2 = [];
        dPhi_dn3 = [];
        dPhi_dnv1x = [];
        dPhi_dnv1y = [];
        dPhi_dnv1z = [];
        dPhi_dnv2x = [];
        dPhi_dnv2y = [];
        dPhi_dnv2z = [];
    end
    
    methods
        function mu = calc_mu(obj, box, mesh, components)                        
            % transform the derivatives of Phi to k-space (these are set in a sub-class)
            obj.dPhi_dn0 = mesh.forward_fft(obj.dPhi_dn0);
            obj.dPhi_dn1 = mesh.forward_fft(obj.dPhi_dn1);
            obj.dPhi_dn2 = mesh.forward_fft(obj.dPhi_dn2);
            obj.dPhi_dn3 = mesh.forward_fft(obj.dPhi_dn3);
            
            obj.dPhi_dnv1x = mesh.forward_fft_vx(obj.dPhi_dnv1x);
            obj.dPhi_dnv1y = mesh.forward_fft_vy(obj.dPhi_dnv1y);
            obj.dPhi_dnv1z = mesh.forward_fft_vz(obj.dPhi_dnv1z);
            
            obj.dPhi_dnv2x = mesh.forward_fft_vx(obj.dPhi_dnv2x);
            obj.dPhi_dnv2y = mesh.forward_fft_vy(obj.dPhi_dnv2y);
            obj.dPhi_dnv2z = mesh.forward_fft_vz(obj.dPhi_dnv2z);
            
            mu = obj.dPhi_dn0.*box.W0 ...
                + obj.dPhi_dn1.*box.W1 ...
                + obj.dPhi_dn2.*box.W2 ...
                + obj.dPhi_dn3.*box.W3 ...
                + obj.dPhi_dnv1x.*-box.Wv1x ... % The vector weighting functions are negative here because of their odd symmetry. See Roth 2010
                + obj.dPhi_dnv1y.*-box.Wv1y ...
                + obj.dPhi_dnv1z.*-box.Wv1z ...
                + obj.dPhi_dnv2x.*-box.Wv2x ...
                + obj.dPhi_dnv2y.*-box.Wv2y ...
                + obj.dPhi_dnv2z.*-box.Wv2z;
            mu = mesh.inverse_fft(mu);
           
        end
        
    end
        
end

