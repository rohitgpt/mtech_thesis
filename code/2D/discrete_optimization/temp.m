% for iter=2:5
nx = 4;
ny = 3;
tic;
E = [4.25e9 20e9];
rho = [1000 2000];
sigma = 5*10e6;
strain = 2*10^(-4);
% force = iter*1e6;

syms xa ya;
%size of each element is such that the overall size is one unit (one centimeter)
%across and one unit wide

alpha = (0.5e-2)/nx;
beta = (0.5e-2)/ny;
n = 1/(4*alpha*beta)*[(alpha-xa)*(beta-ya) (alpha+xa)*(beta-ya) (alpha+xa)*(beta+ya) (alpha-xa)*(beta+ya)];
B = sym('b', [3 8]);
for i=1:4
  B(1:3,i*2-1:i*2) = [diff(n(i), xa) 0; 0 diff(n(i), ya); diff(n(i), ya) diff(n(i),xa)];
end
% B
subs(subs(B, xa, 0), ya, 0)