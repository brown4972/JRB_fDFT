classdef dft_functional_isaft 
    %dft_functional_isaft - iSAFT functional for linear chains
    %   chains must be connected linearly, with each component potentially
    %   bonded to neither, one, or both of its "neighbors", meaning
    %   component j can be bonded to j-1 and/or j+1
    
    
    properties
    end
    
    methods
        function obj = dft_functional_isaft(components)
            if ~isa(components,'dft_polymer_parameters')
                error('ERROR: components must be a dft_polymer_parameters')
            end
        end
        
        function F = calc_free_energy(obj, mu_excess, box, mesh, components)
            % Calculate the free energy from the bonding term only, i.e.
            % equation (8) from Jain et al. 2007, but with the ln X_A terms
            % replaced with the proper propagator terms. the bare X_A
            % term is dropped (it should be small compared to ln X_A, and
            % the 1/2 term goes into the calculation of the free energy 
            % from the chain mu_baths
            [~, P_f, P_b] = calc_fields(obj,mu_excess,box,mesh,components);
            F = sum(box.weighted_rho_seg.*(-P_f-P_b+reshape(components.nGamma/2,1,1,1,components.Nseg)),4);
            F = F.*mesh.cell_weight*mesh.DX^(mesh.D);
            F = sum(F(:));
        end
        
        function [D_field, P_f, P_b] = calc_fields(obj,mu_excess,box,mesh,components)
            [sqrt_y_cav, inv_y_cav, dy_dxi2, dy_dxi3] = obj.calc_y_cav_terms(box, mesh, components);
            
            % NOTE: hard coded linearity here
            y_cav_xi2_term =  inv_y_cav.*dy_dxi2;
            y_cav_xi2_term =  sum(box.weighted_rho_seg.*(y_cav_xi2_term + circshift(y_cav_xi2_term,-1,4)),4);
            
            y_cav_xi3_term =  inv_y_cav.*dy_dxi3;
            y_cav_xi3_term =  sum(box.weighted_rho_seg.*(y_cav_xi3_term + circshift(y_cav_xi3_term,-1,4)),4);
            
            y_cav_xi2_term = box.WXI2.*mesh.forward_fft(y_cav_xi2_term);
            y_cav_xi3_term = box.WXI3.*mesh.forward_fft(y_cav_xi3_term);
            
            D_field = 0.5*mesh.inverse_fft(y_cav_xi2_term + y_cav_xi3_term) - mu_excess;
            
            % NOTE: hard coded linearity here
            P_f = zeros(mesh.NX, mesh.NY, mesh.NZ, components.Nseg);
            P_b = zeros(mesh.NX, mesh.NY, mesh.NZ, components.Nseg);
            seg = 0;
            for poly=1:components.Npoly
                N = components.poly_len(poly);
                seg = seg+1;
                for n=2:N
                    seg = seg+1;
                    
                    % forward walk
                    type_seg_prev = components.seg_to_type(seg-1);
                    integrand = exp(P_f(:,:,:,seg-1)+D_field(:,:,:,type_seg_prev)-box.V_external(:,:,:,type_seg_prev)).*sqrt_y_cav(:,:,:,seg);
                    integrand = mesh.forward_fft(integrand).*box.Wbond(:,:,:,seg);
                    P_f(:,:,:,seg) = real(log(sqrt_y_cav(:,:,:,seg).*mesh.inverse_fft(integrand)));
                    
                    % backward walk
                    seg_b = seg+N-2*n+1;
                    type_seg_prev = components.seg_to_type(seg_b+1);
                    integrand = exp(P_b(:,:,:,seg_b+1)+D_field(:,:,:,type_seg_prev)-box.V_external(:,:,:,type_seg_prev)).*sqrt_y_cav(:,:,:,seg_b+1);
                    integrand = mesh.forward_fft(integrand).*box.Wbond(:,:,:,seg_b+1);
                    P_b(:,:,:,seg_b) = real(log(sqrt_y_cav(:,:,:,seg_b+1).*mesh.inverse_fft(integrand)));
                end
            end
        end
        
        function mu_bath = calc_mu_bath(obj,box,mu_excess,mesh,components)
            % FIXME: This is a bit of a hack, should really do the easier
            % calculaiton of the homogeneous system...
            if ~isequal(box.rho,box.rho0)
                error('ERROR: can only calculate mu_bath if box.rho == box.rho0')
            end
            [D_field, P_f, P_b] = obj.calc_fields(mu_excess,box,mesh,components);
            field_seg = -D_field(:,:,:,components.seg_to_type) -P_f - P_b;
            
            mu_bath = zeros(mesh.NX, mesh.NY, mesh.NZ, components.Npoly);
            for poly = 1:components.Npoly
                mu_bath(:,:,:,poly) = log(components.rho_poly_0(poly)/components.poly_len(poly)) + mean(field_seg(:,:,:,components.seg_to_poly==poly),4);
            end
        end
        
    end
    
    methods (Access = protected)
        % NOTE: hard coded linearity here
        function [sqrt_y_cav, inv_y_cav, dy_dxi2, dy_dxi3] = calc_y_cav_terms(obj, box, mesh, components)
            Delta = 1-box.xi3;
            
            const = zeros(1,1,1,components.Nseg);
            for seg=2:components.Nseg
                if components.seg_to_poly(seg) == components.seg_to_poly(seg-1)
                    d_seg = components.dHS(components.seg_to_type(seg));
                    d_prev = components.dHS(components.seg_to_type(seg-1));
                    const(seg) = d_seg*d_prev/(d_seg+d_prev);
                end
            end
            
            C_xi2_over_Delta = const.*box.xi2./Delta;
            y_cav = ((const>0) + 3*C_xi2_over_Delta + 2*C_xi2_over_Delta.^2)./Delta;
            
            sqrt_y_cav = sqrt(y_cav);
            inv_y_cav = 1./(y_cav);
            inv_y_cav(:,:,:,const==0) = 0;
            
            dy_dxi2 = const.*(3+4*C_xi2_over_Delta)./Delta.^2;
            dy_dxi3 = ((const>0)+6*(C_xi2_over_Delta + C_xi2_over_Delta.^2))./Delta.^2;
        end
    end
    
end

