function [ u ] = u_LJ93_CS( r, m )
%u_LJ93_CS - Lennard-Jones 9-3 potential, cut and shifted at long range
%   m(1) = epsilon
%   m(2) = sigma
%   m(3) = r_cut

u = zeros(size(r));
u_shift = -m(1)*((2/15)*(m(2)/m(3))^9-(m(2)/m(3))^3);
u(r<m(3)) = m(1)*((2/15)*(m(2)./r(r<m(3))).^9-(m(2)./r(r<m(3))).^3) + u_shift;

end

