classdef dft_functional_hard_spheres_asFMT < dft_functional_hard_spheres 
    % dft_functional_hard_spheres_asFMT - the asFMT functional
    %  Fundamental measure theory corrected for 0D crossover, aka "anti-symmetrized"
    %  with numerical stabilization from Knepley 2010
    %
    
    methods
        function F = calc_free_energy(obj, box, mesh, components)
            dot_nv1nv2 = box.nv1x.*box.nv2x + box.nv1y.*box.nv2y + box.nv1z.*box.nv2z;
            dot_nv2nv2 = box.nv2x.^2 + box.nv2y.^2 + box.nv2z.^2;
            dot_nv2nv2_over_n2n2 = min(dot_nv2nv2./(box.n2.^2),1);
            F = -box.n0.*log(1-box.n3) ...
                    + (box.n1.*box.n2 - dot_nv1nv2)./(1-box.n3) ...
                    + (1-dot_nv2nv2_over_n2n2).^3.*box.n2.^3./(24*pi*(1-box.n3).^2);
            F = F.*mesh.cell_weight*mesh.DX^(mesh.D);
            F = sum(F(:));
        end
        
        function mu = calc_mu(obj, box, mesh, components)
            % pre-calculate the dot products and 1-n3
            dot_nv1nv2 = box.nv1x.*box.nv2x + box.nv1y.*box.nv2y + box.nv1z.*box.nv2z;
            dot_nv2nv2 = box.nv2x.^2 + box.nv2y.^2 + box.nv2z.^2;
            dot_nv2nv2_over_n2n2 = min(dot_nv2nv2./(box.n2.^2),1);
            Delta = 1-box.n3;
            
            % calculate the real-space dPhi_dnX's
            obj.dPhi_dn0 = -log(Delta);
            obj.dPhi_dn1 = box.n2./Delta;
            obj.dPhi_dn2 = box.n1./Delta + (1-dot_nv2nv2_over_n2n2).^2.*(box.n2.^2 + dot_nv2nv2)./(8*pi*Delta.^2);
            obj.dPhi_dn3 = box.n0./Delta + (box.n1.*box.n2 - dot_nv1nv2)./(Delta.^2) + (1-dot_nv2nv2_over_n2n2).^3.*box.n2.^3./(12*pi*Delta.^3);
            obj.dPhi_dnv1x = -box.nv2x./Delta;
            obj.dPhi_dnv1y = -box.nv2y./Delta;
            obj.dPhi_dnv1z = -box.nv2z./Delta;
            obj.dPhi_dnv2x = -box.nv1x./Delta - (1-dot_nv2nv2_over_n2n2).^2.*box.n2.*box.nv2x./(4*pi*Delta.^2);
            obj.dPhi_dnv2y = -box.nv1y./Delta - (1-dot_nv2nv2_over_n2n2).^2.*box.n2.*box.nv2y./(4*pi*Delta.^2);
            obj.dPhi_dnv2z = -box.nv1z./Delta - (1-dot_nv2nv2_over_n2n2).^2.*box.n2.*box.nv2z./(4*pi*Delta.^2);
            
            mu = calc_mu@dft_functional_hard_spheres(obj, box, mesh, components);
        end
        
        % uses formula from Knepley et al 2010
        function mu_bath = calc_mu_bath(obj, box, mesh, components)
            dHS = reshape(components.dHS,1,1,1,components.Ncomp);
            
            xi0 = (pi/6)*sum(box.rho0,4);
            xi1 = (pi/6)*sum(box.rho0.*dHS,4);
            xi2 = (pi/6)*sum(box.rho0.*dHS.^2,4);
            xi3 = (pi/6)*sum(box.rho0.*dHS.^3,4);
            Delta = 1-xi3;
            
            P_bath = xi0./Delta + 3*xi1.*xi2./(Delta.^2) + 3*xi2.^3./(Delta.^3);
            
            mu_bath = zeros(size(box.rho0));
            for n=1:components.Ncomp
                mu_bath(:,:,:,n) = -log(Delta) + (3*xi2*dHS(n) + 3*xi1*dHS(n)^2)./Delta ...
                    + 9*xi2.^2*dHS(n)^2./(2*Delta.^2) + P_bath*dHS(n)^3;
            end
        end
    end
        
    
end

