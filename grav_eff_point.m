function [gz] = grav_eff_point(x,xm,m,G)
    G = 6.674*10^(-11);
    r = norm(x-xm);
    gz = G*m*(x(3)-xm(3))/r^3;
end