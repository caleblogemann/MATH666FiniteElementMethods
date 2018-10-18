function [ fh ] = mesh_density( p )
    [theta, rho] = cart2pol(p(:,1), p(:,2));
    a = rho - 1/8;
    b = 0.9 + 0.25*sin(5*theta) - rho;
    c = min(a, b);
    fh = min(c, 2);
end

