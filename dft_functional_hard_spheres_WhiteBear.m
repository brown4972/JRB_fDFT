classdef dft_functional_hard_spheres_WhiteBear < dft_functional_hard_spheres
    % dft_functional_hard_spheres_FMT - the FMT functional
    %  Rosenfeld's original fundamental measure theory 
    %
    
    methods
        function F = calc_free_energy(obj, box, mesh, components)
            dot_nv1nv2 = box.nv1x.*box.nv2x + box.nv1y.*box.nv2y + box.nv1z.*box.nv2z;
            dot_nv2nv2 = box.nv2x.^2 + box.nv2y.^2 + box.nv2z.^2;
            F = -box.n0.*log(1-box.n3) ...
                    + (box.n1.*box.n2 - dot_nv1nv2)./(1-box.n3) ...
                    + (box.n2.^3 - 3*box.n2.*dot_nv2nv2).*(box.n3+(1-box.n3).^2.*log(1-box.n3))./(36*pi*box.n3.^2.*(1-box.n3).^2);
            F = F.*mesh.cell_weight*mesh.DX^(mesh.D);
            F = sum(F(:));
        end
        
        function mu = calc_mu(obj, box, mesh, components)
            % pre-calculate the dot products and 1-n3
            dot_nv1nv2 = box.nv1x.*box.nv2x + box.nv1y.*box.nv2y + box.nv1z.*box.nv2z;
            dot_nv2nv2 = box.nv2x.^2 + box.nv2y.^2 + box.nv2z.^2;
            Delta = 1-box.n3;
            
            % calculate the real-space dPhi_dnX's
            obj.dPhi_dn0 = -log(Delta);
            obj.dPhi_dn1 = box.n2./Delta;
            obj.dPhi_dn2 = box.n1./Delta + (box.n2.^2 - dot_nv2nv2).*(box.n3+Delta.^2.*log(Delta))./(12*pi*box.n3.^2.*Delta.^2);
            obj.dPhi_dn3 = box.n0./Delta + (box.n1.*box.n2 - dot_nv1nv2)./(Delta.^2) ...
                - (box.n2.^3 - 3*box.n2.*dot_nv2nv2).*(box.n3.*(box.n3.^2-5*box.n3+2) + 2*Delta.^3.*log(Delta))./(36*pi*box.n3.^3.*Delta.^3);
            obj.dPhi_dnv1x = -box.nv2x./Delta;
            obj.dPhi_dnv1y = -box.nv2y./Delta;
            obj.dPhi_dnv1z = -box.nv2z./Delta;
            obj.dPhi_dnv2x = -box.nv1x./Delta - box.n2.*box.nv2x.*(box.n3+Delta.^2.*log(Delta))./(6*pi*box.n3.^2.*Delta.^2);
            obj.dPhi_dnv2y = -box.nv1y./Delta - box.n2.*box.nv2y.*(box.n3+Delta.^2.*log(Delta))./(6*pi*box.n3.^2.*Delta.^2);
            obj.dPhi_dnv2z = -box.nv1z./Delta - box.n2.*box.nv2z.*(box.n3+Delta.^2.*log(Delta))./(6*pi*box.n3.^2.*Delta.^2);

            mu = calc_mu@dft_functional_hard_spheres(obj, box, mesh, components);
        end
    end
        
   
   
end

