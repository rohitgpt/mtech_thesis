function main(Problem)
for i=1:7
nx = i;
ny = i;
tic;
E = [20e9 4.25e9];
% rho = E/1e7;
rho = [1300 1000];
% E = [1e9];
% rho = [1000];
p=3;
sigma = 100e6;
strain = 2e-3;
force = 4*1e6/sqrt(i);
[M,min_x, max_x, x] = Problem(nx, ny, E, rho, p, strain, force);
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxFunEvals = 1000000;
options.MaxIter=2000;
options.TolFun=1e-06;
options.ConstraintTolerance = 1e-08;
% disp([1+2*(ny+1)*nx, 2*(ny+1)*(nx+1)-1]);
% disp([3+2*(ny+1)*nx:2:2*(ny+1)*(nx+1)-2]);
[X, ~] = fmincon(@(X) M{1}(X),x,[],[],[],[],min_x,max_x,@(X) mycon(X, M),options);
% M
% disp(1-X);
colormap(gray);
imagesc(1-X, [0, 1]);
title('Material Distribution')
xlabel('X axis')
ylabel('Y axis')
savefig(['same E to rho ratio ', num2str(i),'x',num2str(i)]);
toc
end
end