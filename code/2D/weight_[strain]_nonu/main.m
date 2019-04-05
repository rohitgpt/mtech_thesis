function main(Problem)
nx = 1;
ny = 1;
tic;
% E = [20e9 4.25e9];
% rho = [1300 1000];
E = [1e9];
rho = [1000];
p=1;
sigma = 100;
strain = 2e-3;
[M,min_x, max_x, x] = Problem(nx, ny, E, rho, p, strain);
% options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxFunEvals = 1000000;
options.MaxIter=2000;
options.TolFun=1e-06;
options.ConstraintTolerance = 1e-08;
M
[X, ~] = fmincon(@(X) M{1}(X),x,[],[],[],[],min_x,max_x,@(X) mycon(X, M),options);
disp(1-X);
colormap(gray);
imagesc(1-X, [0, 1]);
% figure(2)
% colormap(gray);
% imagesc(Y(1:a(1), 1:a(2)), [0,1]); 
toc
end