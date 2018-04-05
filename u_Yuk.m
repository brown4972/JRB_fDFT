function [ u ] = u_Yuk( r, m )
%u_Yuk: Yukawa potential
% u(r) = exp(-(r)/L)*E/r 
%   m(1) = E
%   m(2) = L

u = exp(-r/m(2)).*m(1)./r;

end

