function [ u ] = u_1overr4( r, m )
%u_1overr4 - u(r) = -E/(r^4)
%   m(1) = E

u = -m(1)./(r.^4);

end

