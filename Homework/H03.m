f = @(x, y) 8*pi^2*sin(2*pi*x)*sin(2*pi*y);
a = 0;
b = 1;

L2error = [];
Marray = [7, 101, 401];
for M = Marray
    u = uniformSquareCrouziexRaviart(f, M, a, b);
    h = (b-a)/M;
    [X, Y] = meshgrid(a+h/2:h:b-h/2);
    dr = 2*M-1 + M;
    ind = dr*floor(0:1/M:M-1/M) + repmat(1:2:2*M-1,1,M);
    uSol = reshape(u(ind), M, M);
    surf(X, Y, uSol);
    uExact = @(x, y) sin(2*pi*x).*sin(2*pi*y);
    uExactSol = uExact(X, Y);
    L2error = [L2error, norm(uExactSol- uSol)];
end

log(L2error(1:end-1)./L2error(2:end))./log(((b-a)./Marray(1:end-1))./((b-a)./Marray(2:end)))
