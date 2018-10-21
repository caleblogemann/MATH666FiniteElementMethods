%% Problem 1(d)
f = @(x, y) pi^2*exp(sin(pi*x).*sin(pi*y))*(2*sin(pi*x)*sin(pi*y) + 2*cos(pi*x).^2*cos(pi*y).^2 - cos(pi*x).^2 - cos(pi*y).^2);
g = @(x, y) (1 + pi*sin(pi*y))*(x == -1) + (1 - pi*sin(pi*y))*(x==1) + (1 + pi*sin(pi*x))*(y == -1) + (1 - pi*sin(pi*x))*(y == 1);
uExact = @(x, y) exp(sin(pi*x).*sin(pi*y));
N = 20;
a = -1;
b = 1;
u = h02_01(a, b, N, f, g);
uMat = reshape(u, N, N);
h = (b - a)/(N-1);
x = a:h:b;
[X, Y] = meshgrid(x, x);
surf(X, Y,uMat);
view(2);
U = uExact(X, Y);
%surf(X, Y, U);
%view(2);
error = norm(U - uMat)/norm(U)

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

