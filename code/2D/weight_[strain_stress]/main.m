function main(Problem)
for i=2:7
nx = i;
ny = i;
tic;
% E = [20e9 4.25e9];
% rho = E/1e7;
% [1300 1000];
E = [1e9];
rho = [1000];
p=1;
sigma = 20e6;
strain = 2e-3;
force = 1e6/i;
[M,min_x, max_x, x] = Problem(nx, ny, E, rho, p, strain, sigma, force);
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxFunEvals = 1000000;
options.MaxIter=2000;
options.TolFun=1e-06;
% options.ConstraintTolerance = 1e-08;
% disp([1+2*(ny+1)*nx, 2*(ny+1)*(nx+1)-1]);
% disp([3+2*(ny+1)*nx:2:2*(ny+1)*(nx+1)-2]);
[X, ~] = fmincon(@(X) M{1}(X),x,[],[],[],[],min_x,max_x,@(X) mycon(X, M),options);
X
% M{3}(X)
for ely=1:ny
  for elx=1:nx
  double(M{(ely-1)*nx+elx+1}(X))
%   double(M{(ely-1)*2*nx+2*elx+1}(X))
  end
end
% M
% disp(1-X);
colormap(gray);
imagesc(1-X, [0, 1]);
title('Material Distribution')
xlabel('X axis')
ylabel('Y axis')
savefig([num2str(i),'x',num2str(i)]);
toc
end
end