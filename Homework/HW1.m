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
p1 = @(x) -64*pi*(18 - 6*pi^4*x.^2 + 18*pi^2*x.^3 + 2*pi^4*x - 36*pi^2*x.^2 - ...
    27*x + 12*pi^2*x + 4*pi^4*x.^3 + 3*pi^2);
p2 = @(x) 216 - 816*pi^2 - 96*pi^4 + 32*pi^6*x.^2 - 384*pi^4*x.^3 - ...
    64*pi^6*x.^3 - 288*pi^4*x.^2 + 3456*pi^2*x + 576*pi^4*x - 2592*pi^2*x.^2 + ...
    144*pi^4*x.^4 + 32*pi^6*x.^4;
f = @(x) p1(x).*cos(2*pi*x) + p2(x).*sin(2*pi*x);
exactSol = @(x) ((18 + 2*pi^2)*x.^2 + (-24-4*pi^2)*x.^3 + (9 + 2*pi^2)*x.^4).*sin(2*pi*x);
exactSolDD = @(x) 4.*pi.*(2.*(18+2.*pi.^2).*x+3.*((-24)+(-4).*pi.^2).*x.^2+4.*(9+2.* ...
  pi.^2).*x.^3).*cos(2.*pi.*x)+(2.*(18+2.*pi.^2)+6.*((-24)+(-4).* ...
  pi.^2).*x+12.*(9+2.*pi.^2).*x.^2).*sin(2.*pi.*x)+(-4).*pi.^2.*(( ...
  18+2.*pi.^2).*x.^2+((-24)+(-4).*pi.^2).*x.^3+(9+2.*pi.^2).*x.^4).* ...
  sin(2.*pi.*x);
a = 0;
b = 1;

B = @(uDD, vDD) integral(@(x) uDD(x).*vDD(x), a, b, 'AbsTol', 1e-10, 'RelTol', 1e-10);
L = @(v) integral(@(x) f(x).*v(x), a, b,'AbsTol', 1e-10, 'RelTol', 1e-10);

EnergyNorm = @(uDD) sqrt(B(uDD, uDD));
L2Norm = @(u) sqrt(integral(@(x) (u(x)).^2, a, b,'AbsTol', 1e-10, 'RelTol', 1e-10));

hArray = [];
EnergyErrorArray = [];
L2ErrorArray = [];
for M = [9, 19, 39] % 79, 159]
%for M = [9, 19]
    h = (b - a)/(M+1);
    hArray = [hArray, h];
    xj = @(j) a + h*j;
    phi1 = @(x, j) (3*((x-xj(j-1))/h).^2 - 2*((x - xj(j-1))/h).^3).*(x >= xj(j-1) & x < xj(j)) ...
        + (3*((x-xj(j+1))/h).^2 - 2*((x - xj(j+1))/(-1*h)).^3).*(x >= xj(j) & x <= xj(j+1));
    phi2 = @(x, j) ((x - xj(j)).*(x - xj(j-1)).^2/h^2).*(x >= xj(j-1) & x < xj(j)) ...
        + ((x - xj(j)).*(x - xj(j+1)).^2/h^2).*(x >= xj(j) & x <= xj(j+1));
    phi1Mp1 = @(x, j) (3*((x-xj(M))/h).^2 - 2*((x - xj(M))/h).^3).*(x >= xj(M) & x <= xj(M+1));
    phi2Mp1 = @(x) ((x - xj(M+1)).*(x - xj(M)).^2/h^2).*(x >= xj(M) & x <= xj(M+1));
    phi1DD = @(x, j) (6/h^2 - 12*(x - xj(j-1))/h^3).*(x >= xj(j-1) & x < xj(j)) ...
        + (6/h^2 + 12*(x - xj(j+1))/h^3).*(x >= xj(j) & x <= xj(j+1));
    phi2DD = @(x, j) (2*(x - xj(j)) + 4*(x - xj(j-1))).*(x >= xj(j-1) & x < xj(j)) ...
        + (2*(x - xj(j)) + 4*(x - xj(j+1))).*(x >= xj(j) & x <= xj(j+1));
    phi1DDMp1 = @(x) (6/h^2 - 12*(x - xj(M))/h^3).*(x >= xj(M) & x <= xj(M+1));
    phi2DDMp1 = @(x) (2*(x - xj(M+1)) + 4*(x - xj(M))).*(x >= xj(M) & x <= xj(M+1));

    B11ij = @(i, j) integral(@(x) phi1DD(x, i).*phi1DD(x, j) , xj(i-1), xj(i+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B11Mp1j = @(j) integral(@(x) phi1DDMp1(x).*phi1DD(x, j) , xj(M), xj(M+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B11Mp1Mp1 = integral(@(x) phi1DDMp1(x).*phi1DDMp1(x) , xj(M), xj(M+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B12ij = @(i, j) integral(@(x) phi1DD(x, i).*phi2DD(x, j) , xj(i-1), xj(i+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B12Mp1j = @(j) integral(@(x) phi1DDMp1(x).*phi2DD(x, j) , xj(M), xj(M+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B12iMp1 = @(i) integral(@(x) phi1DD(x,i).*phi2DDMp1(x) , xj(M), xj(M+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B12Mp1Mp1 = integral(@(x) phi1DDMp1(x).*phi2DDMp1(x) , xj(M), xj(M+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B22ij = @(i, j) integral(@(x) phi2DD(x, i).*phi2DD(x, j) , xj(i-1), xj(i+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B22Mp1j = @(j) integral(@(x) phi2DDMp1(x).*phi2DD(x, j) , xj(M), xj(M+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    B22Mp1Mp1 = integral(@(x) phi2DDMp1(x).*phi2DDMp1(x) , xj(M), xj(M+1), 'AbsTol', 1e-10, 'RelTol', 1e-10);
    L1i = @(i) integral(@(x) f(x).*(phi1(x, i)), xj(i-1), xj(i+1),'AbsTol', 1e-10, 'RelTol', 1e-10);
    L2i = @(i) integral(@(x) f(x).*(phi2(x, i)), xj(i-1), xj(i+1),'AbsTol', 1e-10, 'RelTol', 1e-10);
    L1Mp1 = integral(@(x) f(x).*(phi1Mp1(x)), xj(M), xj(M+1),'AbsTol', 1e-10, 'RelTol', 1e-10);
    L2Mp1 = integral(@(x) f(x).*(phi2Mp1(x)), xj(M), xj(M+1),'AbsTol', 1e-10, 'RelTol', 1e-10);
    
    xi = cg3(B11ij, B11Mp1j, B11Mp1Mp1, B12ij, B12Mp1j, B12iMp1, B12Mp1Mp1, B22ij, B22Mp1j, B22Mp1Mp1, L1i, L1Mp1, L2i, L2Mp1, a, b, M);

    sol = @(x) sum([cell2mat(arrayfun(@(j) phi1(x, j).*xi(j), 1:M, 'UniformOutput', false)'); phi1Mp1(x)*xi(M+1); cell2mat(arrayfun(@(j) phi2(x, j).*xi(M+1+j), 1:M, 'UniformOutput', false)'); phi2Mp1(x)*xi(end)]);
    solDD = @(x) sum([cell2mat(arrayfun(@(j) phi1DD(x, j).*xi(j), 1:M, 'UniformOutput', false)'); phi1DDMp1(x)*xi(M+1); cell2mat(arrayfun(@(j) phi2DD(x, j).*xi(M+1+j), 1:M, 'UniformOutput', false)'); phi2DDMp1(x)*xi(end)]);
    EnergyError = EnergyNorm(@(x) exactSolDD(x) - solDD(x));
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

    %phi1 = @(x, j) ((1/h^2 - 2*(x - xj(j))/h^3).*(x - xj(j-1)).^2).*(x >= xj(j-1) & x < xj(j)) ...
        %+ ((1/h^2 + 2*(x - xj(j))/h^3).*(x - xj(j+1)).^2).*(x >= xj(j) & x <= xj(j+1));
    %phi2 = @(x, j) ((x - xj(j)).*(x - xj(j-1)).^2/h^2).*(x >= xj(j-1) & x < xj(j)) ...
        %+ ((x - xj(j)).*(x - xj(j+1)).^2/h^2).*(x >= xj(j) & x <= xj(j+1));
    %phi1Mp1 = @(x) ((1/h^2 - 2*(x - xj(M+1))/h^3).*(x - xj(M)).^2).*(x >= xj(M) & x <= xj(M+1));
    %phi2Mp1 = @(x) ((x - xj(M+1)).*(x - xj(M)).^2/h^2).*(x >= xj(M) & x <= xj(M+1));
    %phi1DD = @(x, j) (2*(1/h^2 - 2*(x - xj(j))/h^3) - 8*(x - xj(j-1))/h^3).*(x >= xj(j-1) & x < xj(j)) ...
        %+ (2*(1/h^2 + 2*(x - xj(j))/h^3) + 8*(x - xj(j+1))/h^3).*(x >= xj(j) & x <= xj(j+1));
    %phi2DD = @(x, j) (2*(x - xj(j)) + 4*(x - xj(j-1))).*(x >= xj(j-1) & x < xj(j)) ...
        %+ (2*(x - xj(j)) + 4*(x - xj(j+1))).*(x >= xj(j) & x <= xj(j+1));
    %phi1DDMp1 = @(x) (2*(1/h^2 - 2*(x - xj(M+1))/h^3) - 8*(x - xj(M))/h^3).*(x >= xj(M) & x <= xj(M+1));
    %phi2DDMp1 = @(x) (2*(x - xj(M+1)) + 4*(x - xj(M))).*(x >= xj(M) & x <= xj(M+1));
