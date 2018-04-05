function [ u ] = u_LJ_WCA_0inside( r, m )
%u_LJ_WCA - attractive part of the Lennard-Jones potential, cut and shifted
%long range, and set to 0 inside the hard sphere at short range
%using a Weeks-Chandler-Anderson style approach, for r<r_min, u(r) = u(r_min)
%   m(1) = epsilon
%   m(2) = sigma
%   m(3) = r_cut
%   m(4) = dHS

r_min = 2^(1/6)*m(2);
sigma6 = m(2)^6;
sigma12 = sigma6*sigma6;
u = zeros(size(r));
u_shift = -4*m(1)*(sigma12/m(3)^12 - sigma6/m(3)^6);
u(r<r_min & r>=m(4)) = 4*m(1)*(sigma12/r_min^12 - sigma6/r_min^6) + u_shift;
u(r>=r_min & r<m(3)) = 4*m(1)*(sigma12./r(r>=r_min & r<m(3)).^12 - sigma6./r(r>=r_min & r<m(3)).^6) + u_shift;

end

