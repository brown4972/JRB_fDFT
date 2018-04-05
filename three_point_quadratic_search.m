function x_star = three_point_quadratic_search( f, a, b, fa, fb, min_step )
%find the minimum of f between a and b by evaluating at c=(a+b)/2 and
%fitting quadratically
%   

c = (a+b)/2;
fc = f(c);

% figure(3)
% x = linspace(a,b);
% alpha = -(-a*fc+c*fa+b*fc+a*fb-b*fa-c*fb)/(a*c^2-b*c^2+b*a^2-c*a^2+c*b^2-a*b^2);
% beta = (-a^2*fc+a^2*fb+fa*c^2-fb*c^2+b^2*fc-fa*b^2)/((a-b)*(-a*c+a*b+c^2-b*c));
% zeta = (a^2*b*fc-a^2*c*fb-b^2*a*fc+fb*a*c^2+b^2*c*fa-fa*b*c^2)/((a-b)*(-a*c+a*b+c^2-b*c));
% y = alpha*x.^2 + beta*x + zeta;
% plot([a c b],[fa fc fb],'o',x,y)
% hold on

% From Maple, if we're fitting to f(x)=alpha*x^2+beta*x+zeta
% alpha = (a*(fc-fb)+b*(fa-fc)+c*(fb-fa))/((b-c)*(c-a)*(b-a))
% beta = (a^2*(fb-fc)+b^2*(fc-fa)+c^2*(fa-fb))/((b-c)*(c-a)*(b-a))
% if the parabola is downward facing we may need to cut the domain in half
% and try again, or b may be the optimal choice
% the isnan is also there in case fc turns out to be Inf
alpha_numerator = (a*(fc-fb)+b*(fa-fc)+c*(fb-fa));
while alpha_numerator < 0 || isnan(alpha_numerator)
    % if f(b) is the lowest part of the domain, then it's the best step to
    % take
    if fb<fc && fc<fa
        x_star = b;
%         plot(b,fb,'k*')
%         hold off
%         pause
        return
    end
    
    % otherwise, cut the domain in half and try again
    b = c;
    fb = fc;
    c = (a+b)/2;
    fc = f(c);
%     plot(c,fc,'o')
    alpha_numerator = (a*(fc-fb)+b*(fa-fc)+c*(fb-fa));
end

% the minimum is at the vertex of the parabola
x_star = -(1/2)*(a^2*(fb-fc)+b^2*(fc-fa)+c^2*(fa-fb))/(alpha_numerator);
if x_star < min_step
    x_star = min_step;
elseif x_star > b
    x_star = b;
end

% x = linspace(a,b);
% alpha = -(-a*fc+c*fa+b*fc+a*fb-b*fa-c*fb)/(a*c^2-b*c^2+b*a^2-c*a^2+c*b^2-a*b^2);
% beta = (-a^2*fc+a^2*fb+fa*c^2-fb*c^2+b^2*fc-fa*b^2)/((a-b)*(-a*c+a*b+c^2-b*c));
% zeta = (a^2*b*fc-a^2*c*fb-b^2*a*fc+fb*a*c^2+b^2*c*fa-fa*b*c^2)/((a-b)*(-a*c+a*b+c^2-b*c));
% y = alpha*x.^2 + beta*x + zeta;
% plot([a c b],[fa fc fb],'o',x,y,x_star,f(x_star),'k*')
% hold off
% drawnow

end

