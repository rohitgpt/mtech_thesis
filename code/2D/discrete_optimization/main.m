function main(Problem)
for iter=2:5
nx = 4;
ny = 3;
tic;
E = [4.25e9 20e9];
rho = [1000 2000];
sigma = 5*10e6;
strain = 2*10^(-4);
force = iter*1e6;

syms xa ya;
%size of each element is such that the overall size is one unit (one meter)
%across and one unit wide

alpha = 0.5e-2/nx;
beta = 0.5e-2/ny;
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

[M,min_x, max_x] = Problem(nx, ny, E, Edash, B, rho, strain, sigma, force, KE);
IntCon = 1:nx*ny;
options = optimoptions('ga','UseParallel', true, 'UseVectorized', false);
options.FunctionTolerance = 1e-6;
% options = optimoptions('ga','Algorithm','interior-point');
% options.MaxFunEvals = 1000000;
% options.MaxIter=20000;
% options.TolFun=1e-06;
% options.ConstraintTolerance = 1e-08;
% disp([1+2*(ny+1)*nx, 2*(ny+1)*(nx+1)-1]);
% disp([3+2*(ny+1)*nx:2:2*(ny+1)*(nx+1)-2]);
% M
[X, ~] = ga(@(X) M{1}(X),nx*ny,[],[],[],[],min_x,max_x,@(X) mycon(X, M),IntCon);
X
Xd = zeros(ny, nx);
count=1;
for k = 1:nx
  for j = 1:ny
    Xd(j,k) = X(count);
    count = count + 1;
  end
end
figure(1)
colormap(gray);
imagesc(2-Xd, [0, 2]);
title(['Material Distribution','varying size', num2str(nx),'x',num2str(ny)])
xlabel('X axis')
ylabel('Y axis')
% savefig(['varying max stress ', num2str(sigma/1e6),' Mpa']);
savefig(num2str(iter));
toc
end
end