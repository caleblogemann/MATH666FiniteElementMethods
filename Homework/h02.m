%% Problem 1(d)
f = @(x, y) pi^2*exp(sin(pi*x).*sin(pi*y))*(2*sin(pi*x)*sin(pi*y) + 2*cos(pi*x).^2*cos(pi*y).^2 - cos(pi*x).^2 - cos(pi*y).^2);
g = @(x, y) (1 + pi*sin(pi*y))*(x == -1 && abs(y) < 1) + (1 - pi*sin(pi*y))*(x==1 && abs(y) < 1) + (1 + pi*sin(pi*x))*(y == -1) + (1 - pi*sin(pi*x))*(y == 1);
uExact = @(x, y) exp(sin(pi*x).*sin(pi*y));
a = -1;
b = 1;

hArray = [];
EnergyErrorArray = [];
L2ErrorArray = [];
LInfErrorArray = [];
for N = [10, 20, 40, 80]
    [u, A] = h02_01(a, b, N, f, g);
    uMat = reshape(u, N, N);
    h = (b - a)/(N-1);
    hArray = [hArray, h];
    x = a:h:b;
    [X, Y] = meshgrid(x, x);
    surf(X, Y, uMat);
    view(2);
    U = uExact(X, Y);
    UExact = reshape(U,N^2,1);
    %surf(X, Y, U);
    %view(2);
    e = UExact - u;
    eMat = U - uMat;
    eMatX = [(eMat(2,:) - eMat(1,:))/h; (eMat(3:end, :) - eMat(1:end-2,:))/(2*h); (eMat(end,:) - eMat(end-1,:))/h];
    eMatY = [(eMat(:,2) - eMat(:,1))/h, (eMat(:,3:end) - eMat(:,1:end-2))/(2*h), (eMat(:,end) - eMat(:,end-1))/h];
    EnergyError = norm(h*sqrt(eMatX.^2 + eMatY.^2) + h*eMat, 'fro');
    EnergyErrorArray = [EnergyErrorArray, EnergyError];
    L2Error = norm(h*e);
    L2ErrorArray = [L2ErrorArray, L2Error];
    LInfError = max(abs(e));
    LInfErrorArray = [LInfErrorArray, LInfError];
end
EnergyOrder = log(EnergyErrorArray(1:end-1)./EnergyErrorArray(2:end))./log(hArray(1:end-1)./hArray(2:end));
L2Order = log(L2ErrorArray(2:end)./L2ErrorArray(1:end-1))./log(hArray(2:end)./hArray(1:end-1));
LInfOrder = log(LInfErrorArray(2:end)./LInfErrorArray(1:end-1))./log(hArray(2:end)./hArray(1:end-1));
disp(EnergyOrder);
disp(L2Order);
disp(LInfOrder);

%% Problem 4
h = 0.03;
[p, t, NIN] = sample_mesh(h);
vandermondeMatrix = @(x) [ones(3, 1), x(:,1), x(:,2)];
gradPhiFun = @(phi) (@(j, x, y) [phi(j, 2); phi(j, 3)]);
bilinearForm = @(x, phi1, phi2, gphi1, gphi2, a) a*gphi1(x(1,1), x(1,2))'*gphi2(x(1,1), x(1,2));
f = @(x, y) 1;
u = cg1_2D(p, t, NIN, vandermondeMatrix, gradPhiFun, bilinearForm, f);
trisurf(t, p(:,1), p(:,2), u)
view(2);

