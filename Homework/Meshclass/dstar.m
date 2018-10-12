function d=dstar(p)
    [theta, rho] = cart2pol(p(:,1), p(:,2));
    d=rho - 0.75 - 0.25*sin(5*theta);
end
