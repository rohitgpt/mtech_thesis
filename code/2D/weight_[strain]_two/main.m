function main(Problem)
for i=1:7
nx = i;
ny = i;
tic;
E = [20e9 4.25e9];
rho = [1300 1000];
% E = [1e9];
% rho = [1000];
p=3;
sigma = 50e6;
strain = 2e-4;
force = 2*1e6;

syms xa ya;
alpha = 0.5/nx;
beta = 0.5/ny;
n = 1/(4*alpha*beta)*[(alpha-xa)*(beta-ya) (alpha+xa)*(beta-ya) (alpha+xa)*(beta+ya) (alpha-xa)*(beta+ya)];
B = sym('b', [3 8]);
for i=1:4
  B(1:3,i*2-1:i*2) = [diff(n(i), xa) 0; 0 diff(n(i), ya); diff(n(i), ya) diff(n(i),xa)];
end
nu=0.3;
Edash = 1/(1+nu)/(1-2*nu)*[1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
KE = double(int(int(B'*Edash*B, xa, -alpha, alpha), ya, -beta, beta));
B = subs(subs(B, xa, 0), ya, 0) + ...
subs(subs(B, xa, 1), ya, 0) + ...
subs(subs(B, xa, 0), ya, 1) + ...
subs(subs(B, xa, 1), ya, 1);

[M,min_x, max_x, x] = Problem(nx, ny, E, Edash, B, rho, p, strain, sigma, force, KE);
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxFunEvals = 1000000;
options.MaxIter=20000;
options.TolFun=1e-06;
options.ConstraintTolerance = 1e-08;
% disp([1+2*(ny+1)*nx, 2*(ny+1)*(nx+1)-1]);
% disp([3+2*(ny+1)*nx:2:2*(ny+1)*(nx+1)-2]);
[X, ~] = fmincon(@(X) M{1}(X),x,[],[],[],[],min_x,max_x,@(X) mycon(X, M),options);
% M
% double(M{3}(X))
% disp(1-X);
figure(1)
colormap(gray);
title('Material Distribution')
xlabel('X axis')
ylabel('Y axis')
imagesc(1-X, [0, 1]);
savefig(['mixture|fiber with stress2 ', num2str(nx),'x',num2str(ny)]);
toc
end
end