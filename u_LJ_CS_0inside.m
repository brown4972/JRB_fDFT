function [ u ] = u_LJ_CS_0inside( r, m )
%u_LJ_CS_HS - Lennard-Jones potential, cut and shifted at long range, and
%zero inside the hard sphere at short range
%   m(1) = epsilon
%   m(2) = sigma
%   m(3) = r_cut
%   m(4) = dHS

sigma6 = m(2)^6;
sigma12 = sigma6*sigma6;
u = zeros(size(r));
u_shift = -4*m(1)*(sigma12/m(3)^12 - sigma6/m(3)^6);
u(r>=m(4) & r<m(3)) = 4*m(1)*(sigma12./r(r>=m(4) & r<m(3)).^12 - sigma6./r(r>=m(4) & r<m(3)).^6) + u_shift;

end

