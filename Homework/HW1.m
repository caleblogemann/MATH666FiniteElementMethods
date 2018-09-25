%% Problem 1
a = -pi;
b = pi;
q = @(x) 3 - sin(x) - sin(2*x) - cos(x);
f = @(x) 2*exp(sin(x)).*exp(cos(x));
exactSol = @(x) exp(sin(x)).*exp(cos(x));
exactSolD = @(x) (cos(x) - sin(x)).*exp(sin(x) + cos(x));
B = @(u, uD, v, vD) integral(@(x) uD(x).*vD(x) + q(x).*u(x).*v(x), a, b, 'AbsTol', 1e-10, 'RelTol', 1e-10);
L = @(v) integral(@(x) f(x).*v(x), a, b,'AbsTol', 1e-10, 'RelTol', 1e-10);

EnergyNorm = @(u, uD) sqrt(B(u, uD, u, uD));
L2Norm = @(u) sqrt(integral(@(x) (u(x)).^2, a, b,'AbsTol', 1e-10, 'RelTol', 1e-10));

%     phi = @(x, j) (x - xj(j-1))/(xj(j) - xj(j-1)).*(x >= xj(j-1) & x <= xj(j)) + (x - xj(j+1))/(xj(j) - xj(j+1)).*(x >= xj(j) & x <= xj(j+1));
%     phiD = @(x, j) 1/(xj(j) - xj(j-1)).*(x >= xj(j-1) & x <= xj(j)) + 1/(xj(j) - xj(j+1)).*(x >= xj(j) & x <= xj(j+1));
%     phi0 = @(x) (x - xj(0))/(xj(1) - xj(0)).*(x >= xj(0) & x <= xj(1)) + (x - xj(M+1))/(xj(M) - xj(M+1)).*(x >= xj(M) & x <= xj(M+1));
%     phi0D = @(x) 1/(xj(1) - xj(0)).*(x >= xj(0) & x <= xj(1)) + (x - xj(M+1))/(xj(M) - xj(M+1)).*(x >= xj(M) & x <= xj(M+1));

hArray = [];
EnergyErrorArray = [];
L2ErrorArray = [];
for nCells = [10, 20, 40, 80, 160]
    h = (b - a)/nCells;
    hArray = [hArray, h];
    xj = @(j) a + h*j;
    phi = @(x, j) (x - xj(j-1))/h.*(x >= xj(j-1) & x <= xj(j)) - (x - xj(j+1))/h.*(x >= xj(j) & x <= xj(j+1));
    phiD = @(x, j) 1/h.*(x >= xj(j-1) & x <= xj(j)) - 1/h.*(x >= xj(j) & x <= xj(j+1));
    phi0 = @(x) -(x - xj(1))/h.*(x >= xj(0) & x < xj(1)) + (x - xj(nCells-1))/h.*(x >= xj(nCells-1) & x <= xj(nCells));
    phi0D = @(x) -1/h.*(x >= xj(0) & x < xj(1)) + 1/h.*(x >= xj(nCells-1) & x <= xj(nCells));
    
    Bij = @(i, j) integral(@(x) phiD(x, i).*phiD(x, j) + q(x).*phi(x,i).*phi(x,j), xj(i-1), xj(i+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    Bi0 = @(i) integral(@(x) phiD(x, i).*phi0D(x) + q(x).*phi(x,i).*phi0(x), min(a, xj(i-1)), max(b, xj(i+1)), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B00 = integral(@(x) phi0D(x).*phi0D(x) + q(x).*phi0(x).*phi0(x), xj(0), xj(1), 'AbsTol', 1e-10, 'RelTol', 1e-10) + ...
        integral(@(x) phi0D(x).*phi0D(x) + q(x).*phi0(x).*phi0(x), xj(nCells-1), xj(nCells), 'AbsTol', 1e-10, 'RelTol', 1e-10);

    Li = @(i) integral(@(x) f(x).*(phi(x, i)), xj(i-1), xj(i+1),'AbsTol', 1e-10, 'RelTol', 1e-10);
    L0 = integral(@(x) f(x).*phi0(x), xj(0), xj(1), 'AbsTol', 1e-10, 'RelTol', 1e-10) + ...
        integral(@(x) f(x).*phi0(x), xj(nCells-1), xj(nCells), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    xi = cg1(Bij, Bi0, B00, Li, L0, a, b, nCells);

    sol = @(x) sum([phi0(x)*xi(end); cell2mat(arrayfun(@(j) phi(x, j).*xi(j), 1:nCells-1, 'UniformOutput', false)')]);
    solD = @(x) sum([phi0D(x)*xi(end); cell2mat(arrayfun(@(j) phiD(x, j)*xi(j), 1:nCells-1, 'UniformOutput', false)')]);
    EnergyError = EnergyNorm(@(x) exactSol(x) - sol(x), @(x) exactSolD(x) - solD(x));
    EnergyErrorArray = [EnergyErrorArray, EnergyError];
    L2Error = L2Norm(@(x) exactSol(x) - sol(x));
    L2ErrorArray = [L2ErrorArray, L2Error];
    x = linspace(a, b, nCells*5);
    plot(x, sol(x), x, exactSol(x));
    pause();
end
EnergyOrder = log(EnergyErrorArray(1:end-1)./EnergyErrorArray(2:end))./log(hArray(1:end-1)./hArray(2:end));
L2Order = log(L2ErrorArray(2:end)./L2ErrorArray(1:end-1))./log(hArray(2:end)./hArray(1:end-1));
disp(EnergyOrder);
disp(L2Order);

%% Problem 2
